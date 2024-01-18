# Homework 9
#Paolo Barba, Marta Carriero,Raffaele Liguori.

#We have daily data for the phenomenon of interest that is the rainfall amount, recorded by four different stations (each station has an ID number).
#We have the data from the 1951 to the 1970.
rm(list = ls())
load("datiVE.RData")
str(T4219)
str(T4220)
str(T4237)
str(T4238)



# First of all, since the missing data are a few, in order to don't break the time structure we impute them with the median.
# The missing data are all the same for each of the station.
library(tidyverse)
T4219 <- T4219 %>%
  group_by(year(date)) %>%
  mutate(S4219 = ifelse(is.na(S4219), median(S4219, na.rm = TRUE), S4219))
T4219 <- T4219[,-3]

T4220 <- T4220 %>%
  group_by(year(date)) %>%
  mutate(S4220 = ifelse(is.na(S4220), median(S4220, na.rm = TRUE), S4220))
T4220 <- T4220[,-3]

T4237 <- T4237 %>%
  group_by(year(date)) %>%
  mutate(S4237 = ifelse(is.na(S4237), median(S4237, na.rm = TRUE), S4237))
T4237 <- T4237[,-3]

T4238 <- T4238 %>%
  group_by(year(date)) %>%
  mutate(S4238 = ifelse(is.na(S4238), median(S4238, na.rm = TRUE), S4238))
T4238 <- T4238[,-3]




#We have time series, so we can explore our data graphically in each station.
plot(1:nrow(T4219), T4219$S4219, type = "h", ylab = "Rainfall", main = "T4219")
plot(1:nrow(T4220), T4220$S4220, type = "h", ylab = "Rainfall", main = "T4220")
plot(1:nrow(T4237), T4237$S4237, type = "h", ylab = "Rainfall", main = "T4237")
plot(1:nrow(T4238), T4238$S4238, type = "h", ylab = "Rainfall", main = "T4238")


#1)CHOOSE THE TIME WINDOW TO USE IN MAXIMA BUILDING

# We have to choose the time window in such a way that we obtain iid maxima.
#Since we are dealing with a natural phenomenon, it is a good choice to consider the annaul maxima adn semestral maxima.
# First we try with annual maxima (obtaining iid maxima) and since it produces just 20 maxima per station we move with the semestral time window that produeces 40 maxima that are still iid.

library(lubridate)

max_functions <- function(station){

  max_years_vec <- c()
  for (year in 1951:1970) {
    date_year_1 <- subset(station, year(date) == year & month(date) <= 6)
    max_year_1 <- max(date_year_1[,2], na.rm = TRUE)
    date_year_2 <- subset(station, year(date) == year & month(date) > 6)
    max_year_2 <- max(date_year_2[,2], na.rm = TRUE)
    max_years_vec <- c( max_years_vec , max_year_1 , max_year_2)
  }
  return(max_years_vec)
}
max_functions_year <- function(station){
  max_years_vec <- numeric(20)
  for (year in 1951:1970) {
    date_year <- subset(station, format(date, "%Y") == as.character(year))
    max_year <- max(date_year[,2], na.rm = TRUE)
    max_years_vec[year - 1950] <- max_year
  }
  return(max_years_vec)
}

##FIRST STATION
max_years_19 <- max_functions_year(T4219)
max_semester_19 <- max_functions(T4219)
##SECOND STATION
max_years_20 <- max_functions_year(T4220)
max_semester_20 <- max_functions(T4220)
##THIRD STATION
max_years_37 <- max_functions_year(T4237)
max_semester_37 <- max_functions(T4237)
##FOURTH STATION
max_years_38 <-  max_functions_year(T4238)
max_semester_38 <-  max_functions(T4238)


# Let us see the acf for the maxima considering 1 year time window
par(mfrow = c(2,2))
acf(max_years_19)
acf(max_years_20)
acf(max_years_37)
acf(max_years_38)

# Let us see the acf for the maxima considering 6 months time window
par(mfrow = c(2,2))
acf(max_semester_19)
acf(max_semester_20)
acf(max_semester_37)
acf(max_semester_38)

# Comments: For both one year  and 6 months time-window we obtain uncorrelated maxima for each of the station; only in the second station we can highlight few lags with correlation above the threshold .So we decide to use 6 months as a time-window for our analysis since gives more maxima


plot(max_semester_19, type="o",main= "maxima for T4219 ")
plot(max_semester_20, type="o",main= "maxima for T4220 ")
plot(max_semester_37, type="o",main= "maxima for T4237")
plot(max_semester_38, type="o",main= "maxima for T4237")

##From this plots we observe that we have an uncorrelated identically distributed maxima distribution, that is our goal for the time-window that we choose.


#2)ESTIMATE GEV MODEL FOR EACH STATION


require(ismev)
GEV1 = gev.fit(max_years_19)
#Look at the diagnostic:
gev.diag(GEV1)

GEV2 = gev.fit(max_years_20)
#Look at the diagnostic:
gev.diag(GEV2)

GEV3 = gev.fit(max_years_37)
#Look at the diagnostic:
gev.diag(GEV3)


GEV4 = gev.fit(max_years_38)
#Look at the diagnostic:
gev.diag(GEV4)


#Comments: for the T4219, T4220, T4238 we can see that the probability plot and the quantile plot fit quite well. Also for the return level plot the coverage is satisfied, and the line of the density plot is centered on the empirical histogram.

#For the station T4237 we have different results: the coverage is poor, the quantile plot does not behave so well for large values of the quantiles and the density line is not centered on the empirical histogram.



#COMPUTE RETURN LEVELS OF ORDER 100, 200, 500, 1000:
# N=100, 200,500,1000.

require(evd)
#We can use the function qgev
values= c(1-1/100, 1-1/250, 1-1/500, 1-1/1000)
quantile_19 = qgev(values, loc = GEV1$mle[1],  scale = GEV1$mle[2], shape = GEV1$mle[3])
quantile_20 = qgev(values, loc = GEV2$mle[1],  scale = GEV2$mle[2], shape = GEV2$mle[3])
quantile_37 = qgev(values, loc = GEV3$mle[1],  scale = GEV3$mle[2], shape = GEV3$mle[3])
quantile_38 = qgev(values, loc = GEV4$mle[1],  scale = GEV4$mle[2], shape = GEV4$mle[3])


quantile_19
quantile_20
quantile_37
quantile_38


# Exercise 4
# GPD model:
# on the same data estimate GPD model.
require(ismev)



par(mfrow = c(2,2))
mrl.plot(T4219$S4219, umax = 100)
mrl.plot(T4220$S4220, umax = 100)
mrl.plot(T4237$S4237, umax = 100)
mrl.plot(T4238$S4238, umax = 100)

thr_19 <- 20
thr_20 <- 20
thr_37 <- 25
thr_38 <- 30
par(mfrow = c(1,1))




# GPD model
# The idea behind our study is to set the threshold when the linearity of the Mean Residual Life Plot vanish and we run Maximum Likelihood  estimations with different thresholds (selected visually). Then we will confirm our decision having a look at the diagnostic plots

gpd_19 <- gpd.fit(T4219$S4219 , thr_19)
gpd.diag(gpd_19)
gpd_19$nexc

gpd_20 <- gpd.fit(T4220$S4220 , thr_20)
gpd.diag(gpd_20)
gpd_20$nexc

gpd_37 <- gpd.fit(T4237$S4237 , thr_37)
gpd.diag(gpd_37)
gpd_37$nexc

gpd_38 <- gpd.fit(T4238$S4238 , thr_38)
gpd.diag(gpd_38)
gpd_38$nexc


# Comments: 
# for the fourth station we set an higher threshold since we notice in the density plot a bimodal behavior. In this way we obtain a better result in fitting. 
##Number of exceedances are sufficient to let the estimation procedure for the parameters of the model reliable.


values= c(1-1/ (100 * 365) , 1-1/ (250* 365), 1-1/ (500* 365), 1-1/(1000* 365))
q_19 <- qgpd(values , scale = gpd_19$mle[1],shape = gpd_19$mle[2])
q_20 <- qgpd(values , scale = gpd_20$mle[1],shape = gpd_20$mle[2])
q_37 <- qgpd(values , scale = gpd_37$mle[1],shape = gpd_37$mle[2])
q_38 <- qgpd(values , scale = gpd_38$mle[1],shape = gpd_38$mle[2])

q_19
q_20
q_37
q_38



#6) discuss difference between the two approaches 
# Comparison GEV between and GPD 

##We decide to compare models with BIC because it doesn't be influenced by different sample sizes of models that causes anomalies cases in the likelihood. For GEV model it uses only the maxima that are few points, while for GPD it selects only data above the threshold.
##The likelihood will be higher when n increases, so by using BIC we are penalizing with saample size and estimated parameters.
##For GEV we estimate location, shape and scale.
##For GPD we estimate shape and scale.



QUANTILE_TABLE <- rbind(c(quantile_19,q_19),c(quantile_20,q_20),
             c(quantile_37,q_37),c(quantile_38,q_38))

rownames(QUANTILE_TABLE) <- c("1st return value","2nd return value","3rd return value","4th return value")
colnames(QUANTILE_TABLE) <- c("Quantiles_GEV_19","Quantiles_GPD_19", "Quantiles_GEV_20" , "Quantiles_GPD_20" , "Quantiles_GEV_37", "Quantiles_GPD_37", "Quantiles_GEV_38", "Quantiles_GPD_38")
QUANTILE_TABLE


##For the third station we choose gpd model and it's reasonable with previous diagnostic analysis, since we observe that GEV fits badly with this station.
##Instead for the other stations we have a better result for GEV since their BIC are smaller.