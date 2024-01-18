rm(list=ls())
load("datiAll.RData")
require(spdep)
require(spatialreg)

##Estimate with spautolm and Pseudo-likelihood

#Build the transformed rates variables to regularize them so that they’ll be reasonably modeled using Gaussian distributions. The
#transformation we apply is the square root. 
nc.sids <- nc.sids[nc.sids$CNTY.ID!=2096,]

sids.ft <- sqrt(1000*(nc.sids$SID79/nc.sids$BIR79))+	  
  sqrt(1000*(nc.sids$SID79+1)/nc.sids$BIR79)

#same transformation for the non white born alive
nwbir.ft <- sqrt(1000*(nc.sids$NWBIR79/nc.sids$BIR79))+	   
  sqrt(1000*(nc.sids$NWBIR79+1)/nc.sids$BIR79)		
sids <- data.frame(long=nc.sids$lon,lat=	 nc.sids$lat,sids.ft=sids.ft,nwbir.ft=nwbir.ft,nc.sids)
sids <- sids[sids$CNTY.ID!=2096,]				
#convert the sids data into aspatialobject
coordinates(sids) <- ~long+lat			
coords<-coordinates(sids)				

sids.nb<-dnearneigh(coords, 0, 48.28, longlat=T) #Distance based neighborhood
summary(sids.nb)
sids.W.mat <- nb2mat(sids.nb,style="B", zero.policy=T)                                   
#compute the distances between countries in the neighborhood list
sids.dist <- nbdists(sids.nb,coords,longlat=T)     
#get the vector of distances
dij <- unlist(sids.dist)                           
## elements of the weighting system
term1 <- min(dij)/dij                              
#get id of neighboring counties by row
row.id <- row(sids.W.mat)[sids.W.mat==1]      
# get id of neighboring counties by columns
col.id <- col(sids.W.mat)[sids.W.mat==1]          
#numbers of new borns in 1979
numb <- sids$BIR79                               
## second term of the weigthing system
term2 <- sqrt(numb[row.id]/numb[col.id])          
##Weights
wgts <- term1*term2                               


n <- length(sids$sids.ft)                          # sample size (n=99, Anson has been removed)
#build the weights matrix
wgtmat <- matrix(0,nrow=n,ncol=n)                  

for (i in 1:length(wgts)){
  wgtmat[col.id[i],row.id[i]] <- wgts[i]                            
}
##the matrix must be converted into a weights list to be used in some functions of spdep and related packages 
sids.W <- mat2listw(wgtmat)                        
## this weights system is not symmetrical
isSymmetric(wgtmat)                                

#further more rates depends on the counties birth rates and this induces high heterogeneity
#We need to smooth the variances of the sids rates
# we can modify the weights system that is not symmetrical
#we build a diagonal matrix with the number of new borns per county 

require(stats)
D<-diag(numb,99,99)
wdgtmat<-D%*%wgtmat

isSymmetric(wdgtmat)
sids.WD <- mat2listw(wdgtmat)

sids.WD1 <- nb2listw(sids.nb,style="B", zero.policy=T)                           




###MODELS
##ONLY with intercept  
sids.nullslm <- spautolm(sids.ft ~ 1, data=sids,family="CAR", listw=sids.WD, weights=(BIR79))     
summary(sids.nullslm)
z <- sids.W.mat %*%  sids.ft
pseudo_null <- lm(sids.ft ~ z, data=sids)
summary(pseudo_null)

##### second model: spatial model accounting for the presence of non white new born.  ----
### weights for the least square are given by the total number of children born in 1979
sids.raceslm <- spautolm(sids.ft ~ nwbir.ft, data=sids, family="CAR", listw=sids.WD, weights=(BIR79))                                 
summary(sids.raceslm)                              
pseudo_race <- lm(sids.ft ~ z+nwbir.ft, data=sids, weights = BIR79)
summary(pseudo_race)


LR.Sarlm(sids.nullslm, sids.raceslm)
##We expect that the goodness of fit of the two models differs, as previous analysis in 1974, but comparing them through a likelihood ratio test, data doesn't show enough evidence in order to reject null hypothesis, so two models has the same goodness of fit.
#We can conclude that the covariate "nwbir.ft" is not significant.

##RESIDUAL ANALYSIS
# Assuming you have already fitted models and obtained residuals
residuals_pseudo_null <- residuals(pseudo_null)
residuals_sids_nullslm <- residuals(sids.nullslm)



# Plot residuals against each other
plot(residuals_pseudo_null, residuals_sids_nullslm, 
     main = "Residuals Comparison", 
     xlab = "Pseudo Null Residuals", 
     ylab = "SIDS Nullslm Residuals", 
     col = "darkgreen", pch = 16)
##The residuals show the same behaviour.

plot(sids.ft, type = "l", 
     main = "Real Data VS fitted", 
     xlab = "Index", 
     ylab = "Values", 
     col = "black", lwd = 2)

# Overlay fitted values on the time series plot
lines(fitted(pseudo_null), col = "red", lwd = 2)
lines(fitted(sids.nullslm), col = "blue", lwd = 2)

# Scatter plot of fitted values
plot(fitted(pseudo_null), residuals_pseudo_null, 
     main = "Fitted Values vs. Pseudo Null Residuals", 
     xlab = "Fitted Values", 
     ylab = "Pseudo Null Residuals", 
     col = "red", pch = 16)

# Scatter plot of fitted values
plot(fitted(sids.nullslm), residuals_sids_nullslm, 
     main = "Fitted Values vs. SIDS Nullslm Residuals", 
     xlab = "Fitted Values", 
     ylab = "SIDS Nullslm Residuals", 
     col = "blue", pch = 16)

##Both comparisons between fitted and residuals highlight that they are uncorrelated as we would like to obtain.

#Computing the sum of the squares of the residuals for both models, we can state that the estimates are quite similar.
sum(residuals(pseudo_null)^2) 
sum(residuals(sids.nullslm)^2)

####

sids.raceslm <- spautolm(sids.ft ~ nwbir.ft, data=sids, family="CAR", listw=sids.WD, weights=(BIR79))                                 
pseudo_race <- lm(sids.ft ~ z+nwbir.ft, data=sids, weights = BIR79)

sum(residuals(sids.raceslm)^2) 
sum(residuals(pseudo_race)^2)

##Since we've seen from Likelihood ratio test that the two models have the same goodness of fit, so it's reasonable that residuals of both models show the same results



###CHOOSE THE BEST NEIGHBOROOD 




