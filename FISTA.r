dev.off
rm(list=ls(all=TRUE))
library("glmnet")
setwd("/Users/chenruqian/Desktop/desktop/Academic/WInter\ 2017/STAT538/Syllabus\ and\ Homeworks")

# Input
# Reference for reading inputs http://www.ats.ucla.edu/stat/r/modules/raw_data.htm
myData <- read.table("crime.txt", header = FALSE)

# Pre-processing
# Reference for standardization http://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
df <- myData[,-2]
df[,-1] <- scale(df[,-1])
df[,1] <- df[,1] - mean(df[,1])
# check that we get mean of 0 and sd of 1 for V2, ..., V7 and just mean 0 for V1
colMeans(df) # faster version of apply(scaled.dat, 2, mean)
apply(df, 2, sd)
N <- dim(df)[1] # number of sample points

# Utility functions
L1_norm<- function(x) sum(abs(x)) / length(x)

#################
### Use Package #
###  GLMNET   ###
#################
# Manual https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
# Coercing X into matrix http://stackoverflow.com/questions/30717201/glmnet-not-training-am-i-using-it-wrong
count3 <- 1
lambdaCV <- c()
while(count3 < 5){
  modelCV<- cv.glmnet(as.matrix(df[,-1]), df[,1])
  lambdaCV <- append(lambdaCV, modelCV$lambda.min)
  count3 <- count3 + 1
} 
lambdaCV <- mean(lambdaCV)

#################
betaInfo <- data.frame(matrix(vector(), 1, 5),stringsAsFactors=F)
# fix lambdaCV at the obtained value 34.2453
counter <- 1
seq <- seq(from = 0.5*lambdaCV, to = 1.5*lambdaCV, by = 0.2)
while(counter < 100){
  # bootstrapping, selecting with replacement
  bag <- floor(runif(100, 1, 51)) # each bag has 100 observations
  model<- glmnet(as.matrix(df[bag,-1]), df[bag,1], lambda = rev(seq))
  beta <- as.vector(model$beta[,ceiling(length(seq)/2)])
  betaInfo <- rbind(betaInfo, beta)
  counter <- counter + 1
} 
betaInfo <- betaInfo[-1,]

standardError <- c()
counter2 <- 1
while(counter2 < 6){
  sem<-sd(betaInfo[,counter2])/sqrt(length(betaInfo[,counter2]))
  standardError <- append(standardError,sem)
  counter2 <- counter2+1
}
plot(standardError)



#################
betaInfo <- data.frame(matrix(vector(), 1, 6),stringsAsFactors=F)
# fix lambdaCV at the obtained value 34.2453
counter <- 1
lambdaCVvector <- c()
while(counter < 100){
  # bootstrapping, selecting with replacement
  bag <- floor(runif(100, 1, 51)) # each bag has 100 observations
  modelCVtmp<- cv.glmnet(as.matrix(df[bag,-1]), df[bag,1])
  lambdaCVtmp <- modelCVtmp$lambda.min
  lambdaCVvector <- append(lambdaCVvector, lambdaCVtmp)
  seq <- seq(from = 0.5*lambdaCVtmp, to = 1.5*lambdaCVtmp, by = 0.2)
  model<- glmnet(as.matrix(df[bag,-1]), df[bag,1], lambda = rev(seq))
  beta <- as.vector(model$beta[,ceiling(length(seq))])
  betaInfo <- rbind(betaInfo, beta)
  counter <- counter + 1
} 
betaInfo <- betaInfo[-1,]

standardError <- c()
counter2 <- 1
while(counter2 < 6){
  sem<-sd(betaInfo[,counter2])/sqrt(length(betaInfo[,counter2]))
  standardError <- append(standardError,sem)
  counter2 <- counter2+1
}
plot(standardError)


#############################################################################
dev.off
rm(list=ls(all=TRUE))
library("glmnet")
setwd("/Users/chenruqian/Desktop/desktop/Academic/WInter\ 2017/STAT538/Syllabus\ and\ Homeworks")
# Input
# Reference for reading inputs http://www.ats.ucla.edu/stat/r/modules/raw_data.htm
myData <- read.table("crime.txt", header = FALSE)
# Pre-processing
# Reference for standardization http://stackoverflow.com/questions/15215457/standardize-data-columns-in-r
df <- myData[,-2]
df[,-1] <- scale(df[,-1])
df[,1] <- df[,1] - mean(df[,1])

# Define X, Y matrices 
# X <- as.matrix(myData[,-c(1,2)])
# Y <- as.matrix(myData[,1])
X <- as.matrix(df[,-1])
Y <- as.matrix(df[,1])
N <- dim(X)[1] # number of sample points 
beta <- as.matrix(c(1,2,3,4,5))
lambdaMax <- max(1/50 * abs(t(X)%*%Y)) # = 153.59

softThres <- function(x, lambda) {
  if (x >= lambda) x - lambda
  else if (x<lambda & x > -lambda) 0
  else  lambda + x
}

#vecR <- function(beta, j)  Y-(X[, -j] %*% beta[-j]) # a vector of length N 
#betaUpdate <- function(r, lambda,N,j) {1/N * (X[,j]%*%r), lambda)}

vecR <- function(beta)  Y-(X %*% beta)
betaUpdate <- function(r, lambda,N,j,expre) {if(expre) softThres( (beta[j] + (1/N) * (X[,j]%*%r)), lambda) else beta[j]}

beta <- as.matrix(c(0,0,0,0,0))
lambda <-25

# 154
# 152~100: first coord nonzero
# 50: 1, 3, 4


#vecR1 <- vecR(beta, 1)
#vecR2 <- vecR(beta, 2)
#vecR3 <- vecR(beta, 3)
#vecR4 <- vecR(beta, 4)
#vecR5 <- vecR(beta, 5)

#randomTF <- sample(c(TRUE,FALSE), 1, TRUE)
randInd <- floor(runif(1,1,6))
randomTF <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
randomTF[randInd] <- TRUE

vecRnew <- vecR(beta)

betaUpdate1 <- betaUpdate(vecRnew, lambda,N, 1, randomTF[1])
betaUpdate2 <- betaUpdate(vecRnew, lambda,N, 2, randomTF[2])
betaUpdate3 <- betaUpdate(vecRnew, lambda,N, 3, randomTF[3])
betaUpdate4 <- betaUpdate(vecRnew, lambda,N, 4, randomTF[4])
betaUpdate5 <- betaUpdate(vecRnew, lambda,N, 5, randomTF[5])
beta <- cbind(betaUpdate1,betaUpdate2,betaUpdate3,betaUpdate4,betaUpdate5)

# path wise gradient descent, using previous lambda as warm start

################# scratch ################
softThres2 <- function(x) {
  if (x >= 2) x - 2
  else if (x<2 & x > -2) 0
  else  2 + x
}