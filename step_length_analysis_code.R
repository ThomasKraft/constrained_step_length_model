#Analysis for Venkataraman et al., "Locomotor constraints favor the evolution of the human pygmy phenotype in tropical rainforests"
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(grid)
library(curl)

#download from github ***
url1 <- "https://raw.githubusercontent.com/ThomasKraft/constrained_step_length_model/master/adat2.csv"
adat <- read.csv(curl(url1), as.is=T, header=T)

url2 <- "https://raw.githubusercontent.com/ThomasKraft/constrained_step_length_model/master/wdat2.csv"
wdat <- read.csv(curl(url2), as.is=T, header=T)

# add subject attributes to walking data frame
dat <- merge(wdat, adat, by="id", all.x=T, all.y=F)


# reduce dataset to only relevant columns and walking trials (i.e. remove running), and only to men
dat<-select(dat, id, Condition, Distance, Type, Speed, Time, Num_steps, SL, Population, Sex, Trochanter)

# use only 'N' (Normal speed) as speed. This will result in one value per individual per condition for the Tsimane, and the average of the four trials for the Batek
# remove trail trials, so we only have open and forest
# confine to walking and men
dat<-filter(dat, Speed=="N",Condition!="T",Type=="W",Sex=="M")

# Create velocity column
dat$v<-dat$Distance/dat$Time

## adjust trochanter to meters first
dat$Trochanter<-dat$Trochanter/100

# New column for dimensionless step length
dat$s<-dat$SL/dat$Trochanter # Equation 2

# New column for dimensionless speed
dat$dimV <- dat$v/(9.81*dat$Trochanter)^0.5 # Equation 3

# Make final dataset
final<-dat  %>% group_by(id,Condition,Trochanter,Population) %>% summarise(mSL = mean(SL,na.rm=T), mSL_SD= sd(SL,na.rm=T), mv = mean(v,na.rm=T), mdimV = mean(dimV,na.rm=T), ms = mean(s,na.rm=T))

final<-as.data.frame(final)

# filter out individals who have average dimensionless velocities that are over 0.6 
final <-filter(final, mdimV < 0.6)

#####################################################
##### Calculate speed costs for all individuals #####
#####################################################

cost <- final %>% select(-mSL, -mSL_SD, -mv, -ms) %>% spread(key = Condition, value = mdimV ) %>%
  filter(!is.na(F), !is.na(O)) %>%
  mutate(diff = O-F, speed_cost = diff*(9.8*Trochanter)^0.5)  # back to dimensioned speed

# filter the "final" dataset such that only subjects with both forest and open trials are included

final <- final %>% filter(id %in% cost$id)

final$Condition <- as.factor(final$Condition)
final$Condition <- factor(final$Condition, levels= c("O", "F"))

#### The following are variables that need to be defined early on ###

# Need to define sum_stats here because it's used in Figure 1. This calculates the mean (flat) lines for dimensionless velocity and step length for Figure 4
sum_stats <- final %>% group_by(Population, Condition) %>% summarise(mu.dimV = mean(mdimV), mu.SL = mean(mSL))

## These variables are used for the analysis depicted in figure 5
batek_mean<-sum_stats$mu.dimV[sum_stats$Population == "Batek" & sum_stats$Condition == "O"]
tsimane_mean<-sum_stats$mu.dimV[sum_stats$Population == "Tsimane" & sum_stats$Condition == "O"]

# Batek (mean step length in forest was 0.74)
batek_mSL_forest<- mean(dat$SL[which(dat$Population == "Batek" & dat$Condition == "F")])

# Tsimane (mean step length in forest was 0.6)
tsimane_mSL_forest<- mean(dat$SL[which(dat$Population == "Tsimane" & dat$Condition == "F")])

###################################################################################
#### Figure 1: Theoretical model showing CSLM and anticipate figs 4-5 ####
###################################################################################
# Exponent values for Fig. 1 (NOTE: these values come from empirical models below for Batek, but are used here just to illustrate the general concept. See lines XXX-XXX)
alpha_batek<-exp(0.1356)
beta_batek<-0.2666

fig1a <- ggplot(filter(final,Population=="Batek"),aes(x=Trochanter,y=mdimV,color=Condition)) + 
  xlim(0.7,1) + ylim(0.2,.7) + 
  theme_classic(base_size=18) + 
  ylab(expression(italic(v))) + xlab(expression(italic(L)[leg](m))) + 
  theme(legend.position="none")  + 
  scale_color_manual(values=c( "#E69F00", "green4"))  + 
  geom_segment(aes(x=0.715,xend=0.95,y=sum_stats$mu.dimV[sum_stats$Population == "Batek" & sum_stats$Condition == "O"],yend=sum_stats$mu.dimV[sum_stats$Population == "Batek" & sum_stats$Condition == "O"]),color="#E69F00",lwd=1.5) +
  stat_function(fun = function(x) ((batek_mSL_forest/(alpha_batek*x)))^(1/beta_batek),color="green4",linetype=2, lwd=1.1,xlim=c(0.73,.95)) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  annotate("text",x=0.705,y=.7,label="(A)",size=6)
fig1a

### Next make the speed cost, Fig 1b ###
fig1b<-ggplot(cost[cost$Population=="Batek",],aes(x=Trochanter,y=speed_cost)) +
  xlim(0.73,.95) + ylim(0,.8) +
  theme_classic(base_size=18) +
  xlab(expression(italic(L)[leg](m))) + ylab(expression(paste("Speed cost (", italic(v)[open]," - ",italic(v)[CSLM],")"))) +
  stat_function(fun = function(x) (batek_mean-((batek_mSL_forest/(alpha_batek*x)))^(1/beta_batek))*sqrt(9.81*x),color="grey",linetype=2, lwd=1.2) +  
  annotate("text",x=0.73,y=.8,label="(B)",size=6)+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
fig1b

### Now make figure 1 ###
fig1_final<-grid.arrange(fig1a,fig1b,ncol=2)
ggsave(fig1_final,fileid="fig1_final.tiff",width=12)
# Note that the rest of Figure 1 (arrows and annotations) was made in Powerpoint

############################
### Demographic analysis; see Methods and Materials ###
############################
anthropometry<-adat %>% filter(id %in% unique(cost$id)) %>% group_by(Population) %>% summarise(n = length(unique(id)), mean_weight = mean(Weight, na.rm=T), sd_weight = sd(Weight, na.rm=T), mean_height = mean(Height, na.rm=T), sd_height = sd(Height, na.rm=T), mean_trochanter = mean(Trochanter, na.rm=T), sd_trochanter = sd(Trochanter, na.rm=T), mean_age = mean(Age,na.rm=T), min_LL=min(Trochanter,na.rm=T),max_LL=max(Trochanter, na.rm=T)) 

dat %>% group_by(Population,Sex) %>% summarise(n = length(unique(id)))

anthropometry<-data.frame(anthropometry)

anthropometry

####################################################################################
#### Plot leg length correlation with height ####
####################################################################################
leg_height<-select(adat,Population,Height,Trochanter,Sex)
leg_height<-filter(leg_height,Sex=="M")

mean(filter(leg_height, Population=="Batek")$Trochanter)
# Tsimane
cor.test(leg_height$Trochanter[leg_height$Population=="Tsimane"], leg_height$Height[leg_height$Population=="Tsimane"])
# Pearson's correlation, R = 0.82

# Batek
cor.test(leg_height$Trochanter[leg_height$Population=="Batek"], leg_height$Height[leg_height$Population=="Batek"])
# Pearson's correlation, R = 0.79

###################################################################################
######### Figure 2. Lstep versus Lleg with color coding by environment. Shows that step length is indeed constrained in forest environments ############
###################################################################################
# Batek
fig2a<-ggplot(filter(final,Population=="Batek"),aes(x=Trochanter,y=mSL,color=Condition)) + geom_point(shape=19,size=2.5)  +
  theme_classic(base_size=18) +  xlab(expression(italic(L)[leg](m))) + ylab(expression(italic(L)[step](m))) + theme(legend.position="none")  + scale_color_manual(values=c("#E69F00", "green4")) + stat_smooth(method="lm",size=0.7,se=T, alpha=0.2,aes(fill=Condition)) +
  scale_fill_manual(values=c("#E69F00", "green4")) +
  annotate("text",x=0.805,y=.9,label="(A) Batek",size=6)

# Tsimane
fig2b <- ggplot(filter(final,Population=="Tsimane"),aes(x=Trochanter,y=mSL,color=Condition)) + geom_point(shape=19,size=2.5) + xlim(0.71,.95) + ylim(.5,.85) + theme_classic(base_size=18) +  xlab(expression(italic(L)[leg](m))) + ylab("") + theme(legend.position=c(.225, .765),legend.text=element_text(size=18))  +
  scale_color_manual(values=c("#E69F00", "green4"), labels=c("Open", "Forest"), name="") + stat_smooth(method="lm",size=0.7,se=T, alpha=0.2,aes(fill=Condition)) +
  scale_fill_manual(values=c("#E69F00", "green4"), guide="none") +
  guides(color=guide_legend(override.aes=list(fill=c("#E69F00", "green4")))) +
  annotate("text",x=0.75,y=.845,label="(B) Tsimane",size=6)

##### Make Figure 2 ######
fig2_final<-grid.arrange(fig2a,fig2b,ncol=2)
ggsave(fig2_final,fileid="fig2_final.tiff",width=12)

####### Supporting statistics for figure 2 ###########

#  Batek
summary(lm(mSL~Trochanter,data=final[final$Condition=="O" & final$Population=="Batek",])) # Open regression: F(1,19)=29.55, p < 0.001; Beta = 1.25 +/- 0.23(SE)

summary(lm(mSL~Trochanter,data=final[final$Condition=="F"& final$Population=="Batek",])) # Closed regression: F(1,19)=0.126, p = 0.73; Beta=0.1 +/- 0.28

#  Tsimane
summary(lm(mSL~Trochanter,data=final[final$Condition=="O" & final$Population=="Tsimane",])) # Open regression: F(1,14)=11.57, p = 0.004; Beta=0.85 +/- 0.25(SE)

summary(lm(mSL~Trochanter,data=final[final$Condition=="F"& final$Population=="Tsimane",])) # Closed regression: F(1,14)=0.04, p = 0.85; Beta = -0.04 +/- 0.25

# Do ANCOVA to test for equality of slopes

#Tsimane
anc1_tsi<-aov(mSL~Trochanter*Condition,data=final[final$Population=="Tsimane",])
anc2_tsi<-aov(mSL~Trochanter+Condition,data=final[final$Population=="Tsimane",])

anova(anc1_tsi,anc2_tsi) # Condition is significant for Tsimane (p=0.02)

#Batek
anc1_bat<-aov(mSL~Trochanter*Condition,data=final[final$Population=="Batek",])
anc2_bat<-aov(mSL~Trochanter+Condition,data=final[final$Population=="Batek",])

anova(anc1_bat,anc2_bat) # Condition is significant for Batek (p<0.01)


# T-test on step lengths in open vs forest
# Batek
t.test(final$mSL[final$Condition=="O" & final$Population=="Batek"],final$mSL[final$Condition=="F" & final$Population=="Batek"],paired=T)
# stats: t = 8.02, df = 20, p-value < 0.001, mean of the differences = 0.08

# Tsimane
t.test(final$mSL[final$Condition=="O" & final$Population=="Tsimane"],final$mSL[final$Condition=="F" & final$Population=="Tsimane"],paired=T)
# stats: t = 5.16, df = 15, p-value < 0.001; means difference = 0.09

###################################################################################
######### Figure 3. s versus v. Shows that step lengths follow predicted equation #############################################################################

# Use Equation 1 from text to plot v versus s and get parameter values A and B. Convert to log-log scale to extract parameter values A and B.Regressions based on all data

# For Batek 
summary(lm(log(ms)~log(mdimV),data=filter(final,Population=="Batek")))
# Parameter values: A = e^0.1356 = 1.145, B = 0.2666
# Thus, final equation for Batek is: s = 1.145*v^0.267, R2=0.75

#Calculate residuals
batek_resid<-resid(lm(log(ms)~log(mdimV),data=filter(final,Population=="Batek")))
batek_resid_final<-  data.frame(r=batek_resid,v=log(filter(final,Population=="Batek")$mdimV))
batek_mod<-lm(r~v,data=batek_resid_final)

#Calculate residuals for Tsimane
tsim_resid<-resid(lm(log(ms)~log(mdimV),data=filter(final,Population=="Tsimane")))
tsim_resid_final<-  data.frame(r=tsim_resid,v=log(filter(final,Population=="Tsimane")$mdimV))
tsim_mod<-lm(r~v,data=tsim_resid_final)

# For Tsimane 
summary(lm(log(ms)~log(mdimV),data=filter(final,Population=="Tsimane")))
# Parameter values: A = e^0.16 = 1.17, B = 0.4 
# Thus, final equation for Tsimane is: s = 1.17*v^0.4, R2=0.68

alpha_tsimane<-exp(0.15869)
beta_tsimane<-0.41499

# Batek
fig3a <- ggplot(filter(final,Population=="Batek"),aes(x=mdimV,y=ms,color=Condition)) + geom_point(shape=19,size=2.5)  + xlim(0.2,0.6) + ylim(0.5,1.1) + theme_classic(base_size=18) + ylab(expression(italic(s))) + xlab(expression(italic(v))) + theme(legend.position="none")  + scale_color_manual(values=c("green4", "#E69F00")) +  stat_function(fun = function(x) alpha_batek * (x)^beta_batek,color="darkgrey", lwd=1.25) +  theme_classic(base_size=18) +  theme(legend.position="none")  + scale_color_manual(values=c("#E69F00", "green4")) +
  annotate("text",x=0.27,y=1.05,label="(A) Batek",size=6) +
  geom_segment(aes(x=0.28,xend=0.234,y=0.70,yend=0.77),arrow=arrow(length=unit(0.2,"cm")),color="black") + 
  annotate("text",x=0.3,y=0.685,label="Preferred, eq. 1",size=4.5)


# Tsimane
fig3b <- ggplot(filter(final,Population=="Tsimane"),aes(x=mdimV,y=ms,color=Condition)) +
  geom_point(shape=19,size=2.5) + 
  xlim(0.2,0.6) + ylim(0.5,1.1) +  
  theme_classic(base_size=18) + 
  ylab("") + xlab(expression(italic(v)))  +
  scale_color_manual(values=c("green4", "#E69F00")) +
  stat_function(fun = function(x) alpha_tsimane * (x)^beta_tsimane,color="darkgrey", lwd=1.25) + theme(legend.position=c(.2, .7),legend.text=element_text(size=18)) + scale_color_manual(values=c("#E69F00", "green4"), labels=c("Open", "Forest"), name="")+
  annotate("text",x=0.29,y=1.05,label="(B) Tsimane",size=6) + geom_segment(aes(x=0.24,xend=0.21,y=0.55,yend=0.6),arrow=arrow(length=unit(0.2,"cm")),color="black") + 
  annotate("text",x=0.27,y=0.535,label="Preferred, eq. 1",size=4.5)

##### Make Figure 3 ######
fig3_final<-grid.arrange(fig3a,fig3b,ncol=2)
ggsave(fig3_final,fileid="fig3_final.jpg",width=12)

# t-test on dimensionless velocity between open and forest conditions
# Batek
t.test(final$mdimV[final$Condition=="O" & final$Population=="Batek"],final$mdimV[final$Condition=="F" & final$Population=="Batek"],paired=T)
# stats: t = 18.172, df = 20, p-value < 0.001

# Tsimane
t.test(final$mdimV[final$Condition=="O" & final$Population=="Tsimane"],final$mdimV[final$Condition=="F" & final$Population=="Tsimane"],paired=T)
# stats: t = 8.84, df = 15, p-value < 0.001

##### T-tests on dimensioned velocity in open and closed #####
# Batek
t.test(final$mv[final$Condition=="O" & final$Population=="Batek"],final$mv[final$Condition=="F" & final$Population=="Batek"],paired=T)
# stats: t(20)=17.455, p < 0.001, mean of differences = 0.50 m/s

# Tsimane
t.test(final$mv[final$Condition=="O" & final$Population=="Tsimane"],final$mv[final$Condition=="F" & final$Population=="Tsimane"],paired=T)
# stats: t(15)=8.54, p < 0.001, mean of differences = 0.36 m/s

##############################################################
############################ Figure 4, depicting CSLM ########

fig4a <- ggplot(filter(final,Population=="Batek"),aes(x=Trochanter,y=mdimV,color=Condition)) + xlim(0.7,1) + ylim(0.2,.7) + theme_classic(base_size=18) + ylab(expression(italic(v))) + xlab(expression(italic(L)[leg](m))) + theme(legend.position="none")  + scale_color_manual(values=c( "#E69F00", "green4")) +
  geom_hline(yintercept = sum_stats$mu.dimV[sum_stats$Population == "Batek" & sum_stats$Condition == "O"], colour= "#E69F00", lwd=1.5) +
  stat_function(fun = function(x) ((batek_mSL_forest/(alpha_batek*x)))^(1/beta_batek),color="green4",linetype=2, lwd=1.1,xlim=c(0.74,1)) +
  geom_point(shape=19,size=2.5)+ 
  geom_segment(aes(x=0.95,xend=0.95,y=0.61,yend=0.56),arrow=arrow(length=unit(0.2,"cm")),color="black") + annotate("text",x=0.95,y=0.625,label="Mean v",size=4.5) +
  geom_segment(aes(x=0.95,xend=0.95,y=0.325,yend=0.26),arrow=arrow(length=unit(0.2,"cm")),color="black") +
  annotate("text",x=0.95,y=0.34,label="CSLM",size=4.5)+
  annotate("text",x=0.72,y=.7,label="(A) Batek",size=6)


fig4b <- ggplot(filter(final,Population=="Tsimane"),aes(x=Trochanter,y=mdimV,color=Condition)) + xlim(.7,1)  + ylim(0.2,.7) + theme_classic(base_size=18) + xlab(expression(italic(L)[leg](m))) + ylab("")+
  theme(legend.position=c(.25, .8),legend.text=element_text(size=18)) +
  scale_color_manual(values=c("#E69F00", "green4"), name="", labels=c("Open", "Forest")) +
  geom_hline(yintercept = sum_stats$mu.dimV[sum_stats$Population == "Tsimane" & sum_stats$Condition == "O"], colour= "#E69F00", lwd=1.5) +
  stat_function(fun = function(x) ((tsimane_mSL_forest/x)/alpha_tsimane)^(1/beta_tsimane),color="green4",linetype=2, lwd=1.1)+  geom_point(shape=19, size=2.5) + 
  geom_segment(aes(x=0.96,xend=0.96,y=0.51,yend=0.45),arrow=arrow(length=unit(0.2,"cm")),color="black") + annotate("text",x=0.96,y=0.525,label="Mean v",size=4.5) +
  geom_segment(aes(x=0.96,xend=0.96,y=0.3,yend=0.23),arrow=arrow(length=unit(0.2,"cm")),color="black") +
  annotate("text",x=0.96,y=0.315,label="CSLM",size=4.5) +  
  annotate("text",x=0.73,y=.7,label="(B) Tsimane",size=6)

##### Make Figure 4 ######
fig4_final<-grid.arrange(fig4a,fig4b,ncol=2)
ggsave(fig4_final,fileid="fig4_final.jpg",width=12)

##############################################################
#### Speed cost calculation #############
##############################################################
#### Batek Stats for linear regression####
batek_speed_cost<-lm(speed_cost~Trochanter,data=cost[cost$Population=="Batek",])
summary(batek_speed_cost)
# F(1,19)=12.3, p = 0.002; Equation: y = 2.9388x - 1.96, R2 = 0.36

#### Tsimane stats for linear regression ####
tsimane_speed_cost<-lm(speed_cost~Trochanter,data=cost[cost$Population=="Tsimane",])
summary(tsimane_speed_cost)
# Adjusted R-squared:  0.2211, F-statistic: 5.259 on 1 and 14 DF,  p-value: 0.03783
# Equation is: y = 1.78x - 1.11

### Batek plot for Figure 5A ###
fig5a<-ggplot(cost[cost$Population=="Batek",],aes(x=Trochanter,y=speed_cost)) + xlim(0.73,.95) + ylim(0,.8) + geom_point(shape=19,size=2.5) + theme_classic(base_size=18) + xlab(expression(italic(L)[leg](m))) +   ylab(expression(paste("Speed cost (", italic(v)[open]," - ",italic(v)[CSLM],")")))+ stat_smooth(method="lm",color="black") + stat_function(fun = function(x) (batek_mean-((batek_mSL_forest/(alpha_batek*x)))^(1/beta_batek))*sqrt(9.81*x),color="grey",linetype=2, lwd=1.1) +  
  annotate("text",x=0.75,y=.8,label="(A) Batek",size=6)

### Tsimane plot for Figure 5B ###
fig5b<-ggplot(cost[cost$Population=="Tsimane",],aes(x=Trochanter,y=speed_cost)) + xlim(0.73,.95)  + ylim(0,.8) +geom_point(shape=19,size=2.5) + theme_classic(base_size=18) + xlab(expression(italic(L)[leg](m))) + ylab("")+ stat_smooth(method="lm",color="black") +
  stat_function(fun = function(x) (tsimane_mean-((tsimane_mSL_forest/(alpha_tsimane*x)))^(1/beta_tsimane))*sqrt(9.81*x),color="grey",linetype=2, lwd=1.1)+  
  annotate("text",x=0.755,y=.8,label="(B) Tsimane",size=6)

##### Make Figure 5 ######
fig5_final<-grid.arrange(fig5a,fig5b,ncol=2)
ggsave(fig5_final,fileid="fig5_final.jpg",width=12)

#### RMSE calculations to examine the comparative fit between the linear and the CSLM. See text for details.
RMSE <- function(resids) {
  sqrt(mean((resids)^2))
}

#Batek
#From exponential function
x <- cost$Trochanter[cost$Population == "Batek"]
batek_curve_preds <- (batek_mean-((batek_mSL_forest/(alpha_batek*x)))^(1/beta_batek))*sqrt(9.81*x)
batek_curve_res <- batek_curve_preds - cost$speed_cost[cost$Population == "Batek"]
batek_RMSE_curve <- RMSE(resids=batek_curve_res)

batek_lm_res <- resid(lm(speed_cost~Trochanter,data=cost[cost$Population=="Batek",]))
batek_RMSE_lm <- RMSE(resids = batek_lm_res)

#Tsimane
#From exponential function
x <- cost$Trochanter[cost$Population == "Tsimane"]
tsimane_curve_preds <- (tsimane_mean-((.6/(1.17*x)))^(1/0.4))*sqrt(9.81*x)
tsimane_curve_res <- tsimane_curve_preds - cost$speed_cost[cost$Population == "Tsimane"]
tsimane_RMSE_curve <- RMSE(resids=tsimane_curve_res)

tsimane_lm_res <- resid(lm(speed_cost~Trochanter,data=cost[cost$Population=="Tsimane",]))
tsimane_RMSE_lm <- RMSE(resids = tsimane_lm_res)


######## Now calculate ecological implications of speed cost; see second paragraph of Discussion ########

leg<-c(0.94, 0.84, 0.76) # Suppose leg lengths of industrialized population, Batek, and Efe

dv<-0.5 # Suppose same walking speed in open of dimensionless velocity of 0.5 #

open_v<-dv*sqrt(9.81*leg) ## Walking speeds of  1.518338, 1.435305, and 1.365247 m/s

open_dist<-(open_v*3600*2)/1000 ## Distance covered in two hours at these speeds is pretty similar across statures in the open: 10.93203 10.33419  9.82978 km

### Batek speed cost according to CSLM ###
forest_cost_batek<-(batek_mean-((batek_mSL_forest/(alpha_batek*leg)))^(1/beta_batek))*sqrt(9.81*leg) # Speed cost of 0.91, 0.49, -0.01

forest_speed_batek<-open_v-forest_cost_batek
forest_speed_batek ### 0.5952990, 0.9320699, and 1.3515564 m/s is forest speed for different statures

forest_dist_batek<-(forest_speed_batek*3600*2)/1000 ## Distance covered in two hours at these speeds is drastically different in the forest: 4.286152, 6.710903, and 9.731206 km

### Tsimane speed cost according to CSLM ###
forest_cost_tsimane<-(tsimane_mean-((tsimane_mSL_forest/(alpha_tsimane*leg)))^(1/beta_tsimane))*sqrt(9.81*leg)

forest_speed_tsimane<-open_v-forest_cost_tsimane
forest_speed_tsimane ### 0.8821522, 1.0405820, and 1.2156259 m/s is forest speed for different statures

forest_dist_tsimane<- (forest_speed_tsimane*3600*2)/1000
forest_dist_tsimane ### 6.351496, 7.492190, and 8.752507 km

