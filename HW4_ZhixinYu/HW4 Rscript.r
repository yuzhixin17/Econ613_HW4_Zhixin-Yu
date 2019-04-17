#================================================= R Script for Econ613_HW #4 ==========================================
#================================================= Apirl 10,2019 ==========================================

# import dataset "Koop - Tobias"
ktobias <- read.csv(choose.files(), header = T, na.strings = c("","NA") , stringsAsFactors = F)

#================================================= EXERCISE 1  Data =================================================
#=================================================================
# Question: Represent the panel dimension of wages for 5 randomly selected individuals.
#=================================================================
# first randomly select 5 individuals form 2178 households
id <- sample(1:2178,5)
# substract the first individual's timetrends and wages
id1 <- subset(ktobias[,c(5,3)],ktobias$PERSONID==id[1])
id2 <- subset(ktobias[,c(5,3)],ktobias$PERSONID==id[2])
id3 <- subset(ktobias[,c(5,3)],ktobias$PERSONID==id[3])
id4 <- subset(ktobias[,c(5,3)],ktobias$PERSONID==id[4])
id5 <- subset(ktobias[,c(5,3)],ktobias$PERSONID==id[5])
# we can plot timetrends and wages for these individuals respectively in order to better represent the panel dimension of wages
plot(id1, type = "l")
plot(id2, type = "l")
plot(id3, type = "l")
plot(id4, type = "l")
plot(id5, type = "l")


#================================================= EXERCISE 2  Random Effects =================================================
# Estimate the random effect model under the normality assumption of the disturbance terms
# firstly, substract logwage, educ and potexper as seperate vector
logwage <- as.matrix(ktobias[,3])
educ <- as.matrix(ktobias[,2])
potexper <- as.matrix(ktobias[,4])
# install packages
install.packages("nlme")
library(nlme) 
# fit Linear Model Using Generalized Least Squares, return the estimates
beta_RE <- as.vector(gls(logwage ~ 1 + educ + potexper)$coefficient) 
beta_RE
# [1] 0.79419112 0.09386374 0.03740530


#================================================= EXERCISE 3  Fixed Effects Model =================================================
#=================================================================
# 1/3 Estimate Between estimators
#=================================================================
# calculate mean of logwage, educ and potexper for each individual
Mean <- aggregate(cbind(LOGWAGE,EDUC,POTEXPER) ~ PERSONID, data = ktobias, mean)
# estimate fixed effect
beta_FE_between <- as.vector(lm(Mean[,2] ~ Mean[,3] + Mean[,4])$coefficient)
beta_FE_between
# [1] 0.84556883 0.09309987 0.02599874

#=================================================================
# 2/3 Estimate Within estimators
#=================================================================
# rename colnames and merge data
colnames(Mean) = c("PERSONID","mean_logwage","mean_educ","mean_potexper")
Within <- merge(ktobias, Mean)
# calculate within difference for logwage, educ and potexper
Within$within_logwage <- Within[,3] - Within[,11]
Within$within_educ <- Within[,2] - Within[,12]
Within$within_potexpr <- Within[,4] - Within[,13]
# estimate fixed effect
beta_FE_within <- as.vector(lm(Within[,14] ~ Within[,15] + Within[,16] - 1)$coefficient)
beta_FE_within
# [1] 0.12366202 0.03856107

#=================================================================
# 3/3 Estimate First time difference estimators
#=================================================================
# calculate first difference
Fdiff <- cbind(ktobias[,1],logwage,educ,potexper)
Fdiff1 <- rbind(0,Fdiff)
Fdiff <- rbind(Fdiff,1)
Fdiff2 <- Fdiff-Fdiff1
Fdiff2 <- as.data.frame(Fdiff2[c(1:17920),])
colnames(Fdiff2) <- c("personid","logwage","educ","potexper")
Fdiff2$personid[Fdiff2$personid == 1] = NA
Fdiff2 <- na.omit(Fdiff2)
# estimate fixed effect
beta_FE_Fdiff <- as.vector(lm(Fdiff2[,2] ~ Fdiff2[,3] + Fdiff2[,4])$coefficient)
beta_FE_Fdiff
# [1] 0.049464352 0.038352307 0.003989071


#================================================= EXERCISE 4  Understanding Fixed Effects =================================================
# consider only a random selected 100 individuals
id100 <- sample(1:2178,100)
ktobias100 <- ktobias[ktobias[,1]%in%id100,]

# write probit likelihood function
probit <- function(beta,x,y){
    l = sum(y*log(pnorm(x%*%beta)))-sum((1-y)*log(1-pnorm(x%*%beta)))
    l = -l
    return(l)
}
# substract X and Y vector
X <- as.matrix(ktobias100[,c(2,4)])
Y <- as.matrix(ktobias100[,3])
# estimate the individual fixed effect parameters
beta_IFE <- optim(c(0,0), probit, x=X, y=Y)$par
beta_IFE
# [1] 0.025705738 0.001804506

# run a regression of estimated individual fixed effets on the invariant variables
Mean100 <- aggregate(cbind(LOGWAGE,EDUC,POTEXPER) ~ PERSONID, data = ktobias100, mean)
alpha <- as.matrix(Mean100[,2] - (Mean100[,3]*beta_IFE[1] + Mean100[,4]*beta_IFE[2]))
ktobias100 <- ktobias100[!duplicated(ktobias100$PERSONID),]
X2 <- as.matrix(ktobias100[,6:10])
beta_IFE_iv <- as.vector(lm(alpha ~ X2)$coefficient)
beta_IFE_iv
# [1]  1.836182726 -0.025955048  0.007224095  0.002077072
# [5] -0.243323856  0.003692090

# The standard errors in the previous may not be correctly calculated, because those error terms may correlated across individual fixed effect and are not i.i.d
# In order to robust the result, we should use bootstrap and sample for each individual.


# create an empty matrix boot49
boot49 <- NULL
# write a for loop to calculate each of 49 sample's standard errors
for(i in 1:49){
    # sample rows from datset with replecament
    dat49 <- sample(1:2178,100)
    ktobias49 <- ktobias[ktobias[,1]%in%dat49,]
    # form new X
    X49 <- as.matrix(ktobias49[,c(2,4)])
    # form new Y
    Y49 <- as.matrix(ktobias49[,3])
    boot_beta <- optim(c(0,0), probit, x=X49, y=Y49)$par
    # run regression
    Mean49 <- aggregate(cbind(LOGWAGE,EDUC,POTEXPER) ~ PERSONID, data = ktobias49, mean)
    alpha49 <- as.matrix(Mean49[,2] - (Mean49[,3]*boot_beta[1] + Mean49[,4]*boot_beta[2]))
    ktobias49 <- ktobias49[!duplicated(ktobias49$PERSONID),]
    X49_2 <- as.matrix(ktobias49[,6:10])
    boot_beta_iv <- as.matrix(lm(alpha49 ~ X49_2)$coefficient)
    # save each time result
    boot49 <- rbind(boot49,t(boot_beta_iv))
}    
# calculate stadard errors
boot49_sd <- as.vector(c(sd(boot49[,1]),sd(boot49[,2]),sd(boot49[,3]),sd(boot49[,4]),sd(boot49[,5]),sd(boot49[,6])))
boot49_sd
# [1] 0.15416923 0.05132883 0.01992070 0.01486677 0.12357277 0.01683718
