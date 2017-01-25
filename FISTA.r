##############################################################################################################################
# 
# This file contains ridge regression: gradient descent + accelerated gradient descent (FISTA)
# 
# Implemented FISTA on LASSO with faster convergence rate.
# Idea: instead of just updating (x_t), we now add a new sequence (y_t). Here each y_t is just a linear... 
# ... combination of x_{t-1} and x_t. Then we feed y_t into the update to x_t. 
# Reference paper: A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems∗ Amir Beck† and Marc Teboulle
#
##############################################################################################################################

##########################################################
# Part -2: Clean my workspace and add some libraries 
##########################################################

dev.off
rm(list=ls(all=TRUE))
library("glmnet") # Glmnet is only used here for finding cross-validation choice of lambda. We can do it by hand easily too. 

#####################################
# Part -1: NOTE: set working directory. 
#####################################
setwd("...")

#######################################################################################
# Part 0: Auxiliary - inputs, functions
#######################################################################################

# Input
# Reference for reading inputs http://www.ats.ucla.edu/stat/r/modules/raw_data.htm
# Data source https://web.stanford.edu/~hastie/StatLearnSparsity_files/DATA/crime.html
# Note that this data is not sparse, so we can just do least-squares directly. 
# But even though the data just has 5 predictors, we want to reduce further the number of "important" predictors. 
# Also, the algorithm can be used for sparse models. 

myData <- read.table("crime.txt", header = FALSE) 

X <- as.matrix(myData[,3:7])
Y <- as.matrix(myData[,1])

# starting point
fistaY <- 0
fistaLambda <- 0

# I used glmnet and built-in cross-validation  to find a good lambda.
# Alternatively lambda can be found using cross-validation by hand. I want to upload code by Jan 28th 2017.
modelCV <- cv.glmnet(X, Y, alpha = 0) # alpha=0 represents ridge cost
lambda <- modelCV$lambda.min

# The ridge regression cost function, in penalty form. We want to minimize this over the space of beta. 
f <- function(beta) t(Y- X %*% beta) %*% (Y- X %*% beta) + lambda * t(beta) %*% beta
# The gradient wrt beta is grad_beta = ((Y-X beta)^tr (Y-X beta) + lambda beta^tr beta)  = -2* X^tr (Y-X beta) + 2 * lambda * beta
grad  <- function(beta) -2*t(X) %*% (Y-X %*% beta) + 2*lambda* beta

l2norm<- function(x) sqrt(sum(x^2))

# we can find the least squares beta using Part 1 as well. Just set lambda = 0 and run f, grad again. Then run part 1.
#lambda <- 0
#f <- function(beta) t(Y- X %*% beta) %*% (Y- X %*% beta) + lambda * t(beta) %*% beta
#grad  <- function(beta) -2*t(X) %*% (Y-X %*% beta) + 2*lambda* beta
betaLS <- c(11.2049468, -0.1557923, 13.7888007, 3.6585927, -1.4853808)
betaLSl2norm <- l2norm(betaLS)

#######################################################################################
# Part 1: gradient descent (with armijo to choose step size)
#######################################################################################
beta <- rep(0,5)
t <- 1
valueArray <- c()
betaTraj <- beta
descentStepNum <- 90
while(t<descentStepNum){
  t <- t+1
  # The Armijo backtracking step to find the appropriate step size
  alpha <- 0.5
  gamma <- 0.8
  s <- 1 # step size = 1/L, starting at L=1, i.e. s=1
  delta <- function(beta) (-1) * grad(beta)
  # while loop fixing beta
  fbeta <- f(beta)
  gradel <- t(grad(beta)) %*% delta(beta)
  sArray <- c()
  counter <- 1
  line1 <- c()
  line2 <- c()
  # Below I limited the number of iterations to 100. Even without limiting this, 
  # we experimentally see that generally we need to update s around 60 times.
  # I added an upper bound here for safety, avoid being stuck in a while loop.
  while (f(beta + s*delta(beta)) > fbeta + alpha * s * gradel & counter < 100){
    line1 <- append(line1, f(beta + s*delta(beta)))
    line2 <- append(line2, fbeta + alpha * s * gradel)
    counter <- counter + 1
    sArray <- append(sArray, s)
    s <- gamma * s 
  }
  
  # plot(sArray, main = paste("Armijo tries step size, when t =", t), xlab = "counter", ylab = "s") 
  # In the above plot, we see that the step size decreases as a geometric series (gamma)^n, as defined in the algorithm.
  
  graph <- data.frame(line1, line2)
  matplot(graph,type=c("o","o"), pch=c(1,2),col=c("purple", "red"))
  legend("topright", inset=.05, legend=c("f(beta + s*delta(beta))", "fbeta + alpha * s * gradel"), pch=c(1,2), col=c("purple","red"), horiz=TRUE)
  title(paste("Armijo t =", t))
  
  beta <- beta - s * grad(beta)
  betaTraj <- rbind(betaTraj, t(beta))
  valueArray <- append(valueArray, f(beta))
}
plot(valueArray, main = "ridge cost function", xlab = "step# t", ylab = "cost")
# In the above plot, ideally I see that the cost decreases strictly.

times <- c(0:365)
t1 <- seq(0,365)
plot(times, t1, type="l")


xAxis <- (apply(betaTraj, 1, l2norm) /betaLSl2norm)[-1]
betaGraph <- data.frame(betaTraj[,1], betaTraj[,2],betaTraj[,3],betaTraj[,4],betaTraj[,5]) [-1,]
matplot(xAxis,betaGraph,type=rep("o",5), pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"))
legend("topleft", inset=.05,  legend=c("funding","hs","not-hs","college","college4"),cex = 0.5, pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"), horiz=TRUE)
title(paste("Ridge Regression with num of descents made = ",descentStepNum))



#######################################################################################
# Part 2: FISTA Accelerated gradient descent (with armijo to choose step size)
#######################################################################################

beta <- rep(0,5)
t <- 1
valueArray <- c()
betaTraj <- beta
descentStepNum <- 90
while(t<descentStepNum){
  t <- t+1
  
  alpha <- 0.5
  gamma <- 0.8
  s <- 1 # step size = 1/L, starting at L=1, i.e. s=1
  delta <- function(beta) (-1) * grad(beta)
  # while loop fixing beta
  fbeta <- f(beta)
  gradel <- t(grad(beta)) %*% delta(beta)
  counter <- 1
  sArray <- c()
  while (f(beta + s*delta(beta)) > fbeta + alpha * s * gradel){
    counter <- counter + 1
    sArray <- append(sArray, s)
    s <- gamma * s 
  }
  plot(sArray, main = paste("t =", t), xlab = "counter", ylab = "s") # I want to see the step sizes decreasing
  
  betaNew <- fistaY - s * grad(beta) # beta_(t+1)
  fistaLambdaNew <- 0.5 * (1+ (1+4*fistaLambda^2)^0.5) # lambda_(t+1)
  fistaY <- betaNew + (fistaLambda-1)/fistaLambdaNew * (betaNew - beta) # y_(t+1)
  beta <- betaNew 
  fistaLambda <- fistaLambdaNew
  
  valueArray <- append(valueArray, f(beta))
  betaTraj <- rbind(betaTraj, t(beta))
}
plot(valueArray, main = "accelerated ridge cost function", xlab = "step# t", ylab = "cost")
# In this plot, note that FISTA cost does not decrease strictly, but instead have tiny bumps while decreasing.

xAxis <- (apply(betaTraj, 1, l2norm) /betaLSl2norm)[-1]
betaGraph <- data.frame(betaTraj[,1], betaTraj[,2],betaTraj[,3],betaTraj[,4],betaTraj[,5]) [-1,]
matplot(xAxis,betaGraph,type=rep("o",5), pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"))
legend("bottomright", inset=.05,  legend=c("funding","hs","not-hs","college","college4"),cex = 0.5, pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"), horiz=TRUE)
title(paste("Ridge Regression FISTA with num of descents made = ",descentStepNum))

matplot(betaGraph,type=rep("o",5), main =paste("Ridge Regression FISTA with num of descents made = ",descentStepNum), pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"))
legend("bottomright", inset=.05,  legend=c("funding","hs","not-hs","college","college4"),cex = 0.5, pch=c(1,2,3,4,5),col=c("purple", "red","black","blue","green"), horiz=TRUE)


#######################################################################################
# Part 3: Armijo (backup version)
#######################################################################################
# I separated Armijo backtracking part out in case I need it in the future.
# Armijo algoritm is used to choose step sizes in gradient descent.
# The idea is that we want to go against the gradient to decrease the cost / objective function.
#     The gradient is only a good measurement locally.
#     Therefore we want to take a step that's small so the gradient is still a good measure, but ...
# ... still as large of a step as possible. So we can decrease correctly yet still fast. 
# We do this by starting at a large step size. Then we decrease stepsize gradually (as a geometric sequence at rate "gamma") ...
# ... till we are happy with the decrease (i.e. the condition of the while loop) in the objective function. 
  
alpha <- 0.5
gamma <- 0.8
s <- 1 # step size = 1/L, starting at L=1, i.e. s=1
delta <- function(beta) (-1) * grad(beta)
# while loop fixing beta
fbeta <- f(beta)
gradel <- t(grad(beta)) %*% delta(beta)
counter <- 1
sArray <- c()
while (f(beta + s*delta(beta)) > fbeta + alpha * s * gradel){ # amount of cost function decrease 
  counter <- counter + 1
  sArray <- append(sArray, s)
  s <- gamma * s 
}
plot(sArray, main = paste("t =", t), xlab = "counter", ylab = "s") # I want to see the step sizes decreasing

