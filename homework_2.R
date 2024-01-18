rm(list=ls())
##BAS and CAL are the polygons in metres for bas and cal
load("C:/Users/raffl/Downloads/RainData.RData")
colnames(Basquota32)= colnames(Calabriaquota32)
grid <- rbind(Basquota32,Calabriaquota32)
plot(grid[,c(2,3)], pch=20)
points(dati[dati$anno %in% c(1973,1975),][,c(7,6)], col="red")

rm(list = ls())
load("C:/Users/raffl/Downloads/RainData.RData")
rain_data <- dati[dati$anno %in% c(1973,1975),]
colnames(Basquota32) <- colnames(Calabriaquota32)
grid <- rbind(Basquota32 , Calabriaquota32 )
plot(grid[,c(2,3)],pch = 20)
points(rain_data[,c(7,6)], col = "Red")


# TODO : Explain the structure of the dataset, how many variable and type......


# Variabili quantitative:

summary(rain_data$totanno)

##From the summary we can say this variable has positive asymmetry, since the median is less than the mean and the difference between third quantile and median is higher than the difference between the median and first quantile

# Load the ggplot2 library
library(ggplot2)

# Create a boxplot using ggplot2
ggplot(data = rain_data, aes(x = Nome_Regione, y = totanno)) +
  geom_boxplot(fill = "orange", color = "brown") +
  labs(
    title = "Total Annual Rainfall by Region",
    x = "Region",
    y = "Total Annual Rainfall"
  ) +
  theme_minimal()

##From this graph we notice that Calabria has the highest values of rainfall; actually all the quartiles are higher than the ones of other regions.
##Meaning that in Calabria it rains more often that the others, so this graph leads to us to analyze the reasons why, for example exploring the relationship between rainfall and elevation



# Create a boxplot using ggplot2
ggplot(data = rain_data, aes(x = Nome_Regione, y = Quota)) +
  geom_boxplot(fill = "orange", color = "brown") +
  labs(
    title = "Elevation by Region",
    x = "Region",
    y = "Elevation"
  ) +
  theme_minimal()

##We observe the elevation for the 3 different regions and we notice that Calabria has a less value for 1st quartile and median respect than Basilicata, but 3nd quartile is higher and it has an heavy code, basically we can say that it has more variability.


# Create a scatter plot using ggplot2 with axis labels, a title, and reference lines
ggplot(data = rain_data, aes(x = Quota, y = totanno)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Scatter Plot of Quota vs. Total Annual Rainfall",
    x = "Quota (Elevation)",
    y = "Total Annual Rainfall (mm)"
  ) +
  theme_minimal()
cor(rain_data$totanno, rain_data$Quota)
##Now we go to analyze the relationship between Elevation and Rainfall with a Scatterplot. We notice that there isn't an high correlation, in fact correlation is 0.39


ggplot(data = rain_data, aes(x = Quota, y = totanno, color = Nome_Regione)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Scatter Plot of Quota vs. Total Annual Rainfall",
    x = "Quota (Elevation)",
    y = "Total Annual Rainfall (mm)"
  ) +
  theme_minimal() +
  facet_wrap(~Nome_Regione)

cor(rain_data$totanno[rain_data$Nome_Regione=="Basilicata"],rain_data$Quota[rain_data$Nome_Regione=="Basilicata"])
cor(rain_data$totanno[rain_data$Nome_Regione=="Calabria"],rain_data$Quota[rain_data$Nome_Regione=="Calabria"])
cor(rain_data$totanno[rain_data$Nome_Regione=="Puglia"],rain_data$Quota[rain_data$Nome_Regione=="Puglia"])

##Commentare questi grafici



# VAriabili qualitative

# Boxplot di totanno condizionato alla regione


# Scatter plot tra quota totatanno

# 3 scatter plot divisi per regione 




# Es 2 ------------------------------------------------------------------------------------------------------------------------------------------

# Riusare la funzione per cross validation
idx <- sample(1:nrow(rain_data), 0.7 * nrow(rain_data))  # Generate a random sample of indexes
train <- rain_data[idx,]   # Make the train set
test <- rain_data[-idx,]   # Make the test set



cross_val_func <- function(data, k , p){
  
  coordinates(data) <- ~Lon10+Lat10  # Transform the data into coordinates format
  
  l_folds <- floor(length(data) / k)  # Compute the length of the fold
  mse <- numeric(k)   # Pre set the mse vector
  
  # Loop over the fold
  for (i in 1:k) {
    ind <- ((i-1)*l_folds+1): (i*l_folds)   # compute the needed index
    cv_test <- data[ind,]    # Training data
    cv_train <- data[-ind,]  # Validation data
    
    # train the model and compute the predictions
    predictions <- idw(formula = totanno ~ 1, locations = cv_train , newdata = cv_test,  idp = p)$var1.pred
    # Compute the mse
    mse[i] <- mean((predictions - cv_test$totanno)^2) 
    
  }
  
  d <- c(mean(mse) , sd(mse))  # Store the mean and sd of the mse's
  return(d) # return mean and sd of the mse of the 5 folds
}


choose_p <- function(grid_p,data , k){
  mq <- c()   # Preset the our proper score vector
  # Loop over the grid parameters
  for (x in grid_p) {
    a <- cross_val_func(data,k,x) # Perform  cross validation 
    mq <- c( mq, a[1] +a[2])      # Store the results
  }
  return(c(grid_p[which.min(mq)] ,  min(mq) ))
  
}
par <- seq(0,10,by = 0.1) # Grid of parameter used

# Perform  5-fold cross validation
best_p <- choose_p(grid_p = par, data = train, k = 5) 
p_opt <- best_p[1]  # Select the best p
p_opt





# ES -----------------------------------------------------------------------------------------------------------------------------------------------


# Usando le stime di rainfall fare un bel grafico
