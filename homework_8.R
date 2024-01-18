#HOMEWORK 8
#Using the data in the mite datasets (from library vegan in R) create the presence-absence variable for species ONOV and LCIL
rm(list = ls())
library(vegan)
library(sp)
library(spdep)
data(mite)
data(mite.env)
data(mite.pcnm)
data(mite.xy)

dat <- mite.env
#PRESENE_ABSENCE
ifelse(mite$ONOV!=0,1,0)
dat$ONOV=ifelse(mite$ONOV!=0,1,0)
dat$LCIL=ifelse(mite$LCIL!=0,1,0)
#It is a dataframe that contains 70 data points and 35 variables that are species of mites. We are interested in 2 specific species that are ONOV and LCIL



#We set the coordinates 
dat$x= mite.xy$x
dat$y=mite.xy$y
coordinates(dat)= ~x+y




#1. Explore the available variables (environmental variables are in the mite.env)
str(mite.env)
#The environmental variables contained in mite.env are:
#-SubsDens:  Substrate density
#-WatrCont:  Water content of the substrate
#-Substrate: Substrate type 
#-Topo:      Microtopography
#The first 2 variables are numeric, the last 3 are categorial variables.
chisq.test(dat$ONOV, dat$Topo, simulate.p.value = T) 
chisq.test(dat$ONOV, dat$Substrate, simulate.p.value = T) 
chisq.test(dat$ONOV, dat$Shrub, simulate.p.value = T) 
###so we're going to consider shrub and substrate as covariates of our autologistic model since data show enough evidence in order to  reject the null hypothesis of independence between the two variables and to not reject it for "Topo" variable.

t.test(dat$SubsDens~dat$ONOV) 
t.test(dat$WatrCont~dat$ONOV) 
t.test(dat$SubsDens~dat$LCIL) 
t.test(dat$WatrCont~dat$LCIL) 

##We're going to consider also this two quantitative variables since we reject null hypothesis for ONOV.
##For LCIL we don't reject the null hypothesis for "SubsDens" variable, so we decide to not include it in the model.
#For "WatrCont" we reject null hypothesis so we will consider it in our models. 


##ANALYSIS FOR ONOV

# Exercise 2: choose the best neighboord structure

coords <- coordinates(dat)


# Build the autogistic model
distances <- seq(20,80 , by = 5)
aic <- rep(NA , length(distances))
i <- 0
for (d in distances){
  i <- i+1
  sids.nb<-dnearneigh(coords, 0, d, longlat=T)
  sids.W.mat <- nb2mat(sids.nb, style="B", zero.policy=T)
  Z <- sids.W.mat %*% dat$ONOV
  
  
  autologistic_model <- glm(ONOV ~  Z  + SubsDens + WatrCont  -1 ,
                            data = dat,
                            family = binomial(link = "logit"))
  aic[i] <- autologistic_model$aic
  
}

d <- distances[which.min(aic)]
sids.nb<-dnearneigh(coords, 0, d, longlat=T)
sids.W.mat <- nb2mat(sids.nb, style="B", zero.policy=T)
Z <- sids.W.mat %*% dat$ONOV


autologistic_complete <- glm(ONOV ~  Z  + SubsDens + WatrCont  -1 ,
                          data = dat,
                          family = binomial(link = "logit"))

##We choose the best distance for the autologistic models with all quantitative covariates and we select the distance with the lowest AIC which correspond to 65
##For simplicity we use only the best distance computed for the complete model.

plot(aic ,
     type = "o", 
     lwd = 2, 
     col = "red",
     main = "AIC with model considering \n different neighboord strctures" , 
     xlab = "Distance", axes = FALSE)
axis(1, 1:length(aic), labels = distances)  #
axis(2)




##Before drawing the predictive map, we compute the best threshold for different models in order to maximize the balanced score, since it is a unbalanced dataset
library(caret)
get_the_best_threshold <- function(model , y){
  thresholds <- seq(.55,.75,by = .05)
  scores <- rep(NA , length(thresholds))
  for (t in 1:length(thresholds)) {
    pp <- predict(model ,type="response")
    fit <- ifelse(pp> thresholds[t],1,0)
    scores[t] <-  confusionMatrix(factor(y) ,factor(fit, levels=c(0,1)))[["byClass"]][11] 
    
  }
  best_score <- max(scores, na.rm = T)
  best_threshold <- thresholds[which.max(scores)]
  return(list(Best_score = best_score , Best_threshold =   best_threshold ))
}
autologistic_simple= glm(ONOV ~  Z  -1 ,
                         data = dat,
                         family = binomial(link = "logit"))

autologistic_one_covariate <- glm(ONOV~Z-1+WatrCont, data=dat, family = binomial(link="logit"))

threshold_simple <- get_the_best_threshold(autologistic_simple , dat$ONOV)
threshold_complete<-  get_the_best_threshold(autologistic_complete, dat$ONOV)
threshold_one <- get_the_best_threshold(autologistic_one_covariate,dat$ONOV)


##PREDICTIVE MAPS

##SIMPLE MODEL
pred.nspat<-fitted(autologistic_simple)
pred.nspat_simple<-ifelse (pred.nspat>=threshold_simple$Best_threshold,1,0)
x <- mite.xy$x
y <- mite.xy$y
plot(x,y,type="n", xlab="X", ylab="Y")
text(x, y, format(round(Z, 0)), pos = 2) 
for (i in 1:70) {
  if (pred.nspat_simple[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}





title("ONOV data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)


##ONE COVARIATE MODEL

pred.nspat<-fitted(autologistic_one_covariate)
pred.nspat_one<-ifelse (pred.nspat>=threshold_one$Best_threshold,1,0)

plot(x,y,type="n", xlab="X", ylab="Y")
text(x,y, format(round(Z,0)), pos=2)

for (i in 1:70) {
  if (pred.nspat_one[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}




title("ONOV data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)

###COMPLETE MODEL
pred.nspat<-fitted(autologistic_complete)
pred.nspat_complete<-ifelse (pred.nspat>=threshold_complete$Best_threshold,1,0)

plot(x,y,type="n", xlab="X", ylab="Y")
text(x,y, format(round(Z,0)), pos=2)

for (i in 1:70) {
  if (pred.nspat_complete[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}




title("ONOV data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)


###BEST PREDICTIVE MODELS CHOICE
cbind(threshold_complete$Best_score, threshold_one$Best_score, threshold_simple$Best_score)
##We have similar results in predictive ability for simple and one covariate model, while the best is the complete model.




####ANALYSIS FOR LCIV variable

Z <- sids.W.mat %*% dat$LCIL

# Build the autologistic model
distances <- seq(20,80 , by = 5)
aic <- rep(NA , length(distances))
i <- 0
for (d in distances){
  i <- i+1
  sids.nb<-dnearneigh(coords, 0, d, longlat=T)
  sids.W.mat <- nb2mat(sids.nb, style="B", zero.policy=T)
  Z <- sids.W.mat %*% dat$LCIL
  
  
  autologistic_model <- glm(LCIL ~  Z  + SubsDens + WatrCont  -1 ,
                            data = dat,
                            family = binomial(link = "logit"))
  aic[i] <- autologistic_model$aic
  
}
d <- distances[which.min(aic)]
d
##We choose the best distance for the autologistic models with all quantitative covariates and we select the distance with the lowest AIC which correspond to 60
##For simplicity we use only the best distance computed for the complete model.

sids.nb<-dnearneigh(coords, 0, d, longlat=T)
sids.W.mat <- nb2mat(sids.nb, style="B", zero.policy=T)
Z <- sids.W.mat %*% dat$LCIL


autologistic_complete <- glm(LCIL ~  Z  + SubsDens + WatrCont  -1 ,
                             data = dat,
                             family = binomial(link = "logit"))


plot(aic ,
     type = "o", 
     lwd = 2, 
     col = "red",
     main = "AIC with model considering \n different neighboord strctures" , 
     xlab = "Distance", axes = FALSE)
axis(1, 1:length(aic), labels = distances)  #
axis(2)




##As we did  before, we compute the best threshold for different models in order to maximize the balanced score, since it is a unbalanced dataset
autologistic_simple= glm(LCIL ~  Z  -1 ,
                         data = dat,
                         family = binomial(link = "logit"))

autologistic_one_covariate <- glm(LCIL~Z-1+WatrCont, data=dat, family = binomial(link="logit"))

threshold_simple <- get_the_best_threshold(autologistic_simple , dat$LCIL)
threshold_complete<-  get_the_best_threshold(autologistic_complete, dat$LCIL)
threshold_one <- get_the_best_threshold(autologistic_one_covariate,dat$LCIL)


##PREDICTIVE MAPS 
##SIMPLE MODEL
pred.nspat<-fitted(autologistic_simple)
pred.nspat_simple<-ifelse (pred.nspat>=threshold_simple$Best_threshold,1,0)

plot(x,y,type="n", xlab="X", ylab="Y")
text(x,y, format(round(Z,0)), pos=2)

for (i in 1:70) {
  if (pred.nspat_simple[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}




title("LCIL data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)


##ONE COVARIATE MODEL

pred.nspat<-fitted(autologistic_one_covariate)
pred.nspat_one<-ifelse (pred.nspat>=threshold_one$Best_threshold,1,0)

plot(x,y,type="n", xlab="X", ylab="Y")
text(x,y, format(round(Z,0)), pos=2)

for (i in 1:70) {
  if (pred.nspat_one[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}




title("LCIL data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)

###COMPLETE MODEL
pred.nspat<-fitted(autologistic_complete)
pred.nspat_complete<-ifelse (pred.nspat>=threshold_complete$Best_threshold,1,0)

plot(x,y,type="n", xlab="X", ylab="Y")
text(x,y, format(round(Z,0)), pos=2)

for (i in 1:70) {
  if (pred.nspat_complete[i] == 1) 
    points(x[i], y[i], col = "red", pch =16, cex = 2, lwd=0.2)
}




title("LCIL data posting", cex.main=1.5)
mtext(side=3, line=0.3, "Circles highlight estimated presence",
      cex=1.0)


###BEST PREDICTIVE MODELS CHOICE
cbind(threshold_complete$Best_score, threshold_one$Best_score, threshold_simple$Best_score)

##As we expected the model with only the "WatrCont" variable is better since in data exploratory analysis we show that "SubsDens" was no informative for LCIV.
