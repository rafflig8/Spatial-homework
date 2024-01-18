################################
################################
########## Homework 1 ########## 
################################
################################

#  Authors : Marta  Carriero , Raffaele Liguori , Paolo Barba
#  Course  : Spatial Statistics and Statistical Tools for Environmental Data
rm(list = ls())
set.seed(123)
load("datiAll.RData")
packages <- c("geoR" ,  "gstat" , "sp" , "MBA" , "viridis")
invisible(lapply(packages, require, character.only = TRUE))


# Es 1 ----------------------------------------------------------------------------------------------------------


idx <- sample(1:nrow(wolf), 0.7 * nrow(wolf))  # Generate a random sample of indexes
train <- wolf[idx,]   # Make the train set
test <- wolf[-idx,]   # Make the test set

# Function to perform k-fold cross validation from a given data set.
# The model selected is the inverse  weighting distance one.
cross_val_func <- function(data, k , p){
  
  coordinates(data) <- ~x+y  # Transform the data into coordinates format
  
  l_folds <- floor(length(data) / k)  # Compute the length of the fold
  mse <- numeric(k)   # Pre set the mse vector
  
  # Loop over the fold
  for (i in 1:k) {
    ind <- ((i-1)*l_folds+1): (i*l_folds)   # compute the needed index
    cv_test <- data[ind,]    # Training data
    cv_train <- data[-ind,]  # Validation data
    
    # train the model and compute the predictions
    predictions <- idw(formula = piezo ~ 1, locations = cv_train , newdata = cv_test,  idp = p)$var1.pred
    # Compute the mse
    mse[i] <- mean((predictions - cv_test$piezo)^2) 

  }
  
  d <- c(mean(mse) , sd(mse))  # Store the mean and sd of the mse's
  return(d) # return mean and sd of the mse of the 5 folds
}

# Below we define the "choose_p" function.
# It is a function in order  to choose the best p parameter of the inverse weighting distance model
# We define our proper score as the sum between fold's mean mse and the fold's sd for each of the parameter selected. We will choose as the best parameter the parameter associated with the minimum value of that score.
# We select this score since with want to control both the standard deviation and the mean of our estimation.

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


# Es 2 ----------------------------------------------------------------------------------------------------------


# Create the gird 20X20
x <- seq(floor(min(train$x)), floor(max(train$x)) , length.out = 20)
y <- seq(floor(min(train$y)), floor(max(train$y)) , length.out = 20)

d1 <- expand.grid(x = x, y = y) # Expand the grid
coordinates(train) <- ~x+y      # Proper set-up
gridded(d1) <- ~x+y             # Grid

# Train the model and perform over the grid
grid_estimation <- idw(formula = piezo ~ 1, locations = train , newdata = d1,  idp = p_opt)

my_palette <- colorRampPalette(c("blue", "yellow", "red"))(100) # Choose the palette

# Plot
spplot(
  oo["var1.pred"],
  col.regions = my_palette, 
  main = list(label = "Predicted Values of piezo \n from x and y coordinates"), 
  colorkey = list(space = "left", height = 0.4), 
  scales = list(draw = TRUE), 
  col = "transparent"
)

# Comments about the plot: 
# From the graphical representation we can see that there are basically 3 clusters  well defined. In the bottom-left there are the observations with the highest values of the target "piezo" variable. Going up in the graph (around the diagonal) we notice  decreasing values of the variable. On the top right we can remark the third class, which present the lowest value of the "piezo" variable.


# Here we want to compare two selected model:
# idw and MBA ones. The comparing measure we will use is the mse.


coordinates(test) <- ~x+y  # Set-up the coordinates

# Train the idw model and predict with the test observations
predictions <- idw(formula = piezo ~ 1, locations = train , newdata = test,  idp = p_opt)$var1.pred
# Compute the mse for idw
mse_idw <- mean((predictions - test$piezo)^2)
# Train the mba model and predict with the test observations
xyz <- data.frame(as.numeric(train$x) , as.numeric(train$y) , as.numeric(train$piezo))
pred_mba <- mba.points(xyz , coordinates(test))$xyz.est[,3]
# The output of this command line returns us a warning message since it can be possible that one or more observations inside the test set is outside the range imposed by the train grid (defined with the extreme value of observations).
# Compute the mse for the mba model
mse_mba <- mean((pred_mba - test$piezo)^2)


# Compare the two mse
mse_idw
mse_mba

# Results: From our case study we conclude that the model that best interpolate our data on a 20x20 grid is the mba.




