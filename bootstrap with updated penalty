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
#either eventually get lambdaCV 0.1706233
df <- scale(df)
#or eventually get lambdaCV 50.1528
#df[,-1] <- scale(df[,-1])
#df[,1] <- df[,1]-mean(df[,1])

# check that we get mean of 0 and sd of 1 for V2, ..., V7 and just mean 0 for V1
colMeans(df) # faster version of apply(scaled.dat, 2, mean)
apply(df, 2, sd)
N <- dim(df)[1] # number of sample points
X <- as.matrix(df[,-1])
Y <- as.matrix(df[,1])

# Utility functions
l1norm<- function(x) sum(abs(x))
r <- function(x, y, beta) y - x%*% beta 
softThres <- function(x, lambda) {
  if (x >= lambda) x - lambda
  else if (x<lambda & x > -lambda) 0
  else  lambda + x
}
lassoCost <- function(beta,x,y,lambda) {t(y- x %*% beta) %*% (y- x %*% beta) + lambda *l1norm(beta) }


#######################################
### cross validation to find lambda
#######################################

bag <- floor(runif(bootstrapBagSize, 1, 51)) # probability that the bag has entires=51 is zero
beta <- rep(0,5) # c(0,0,0,0,0)
perm <- sample(c(1:50), 50, replace = TRUE, prob = NULL)
X <- X[perm,]
Y <- Y[perm]
Xtrain <- X[c(1:40),]
Ytrain <- Y[c(1:40)]
N <- dim(Xtrain)[1]
Xtest <- X[c(41:50),]
Ytest <- Y[c(41:50)]
lengthLambda <- 50
lambdaMax <- max(1/N* t(Xtrain)%*%Ytrain)
lambdaVec <- rev(seq(from = 0, to = lambdaMax, by = (lambdaMax -0)/(lengthLambda-1)))


beta <- rep(0,5)
lassoCostLambdaPath <- c()

for (j in c(1:lengthLambda)) {
  
  betaTraj <- beta
  oldNorm <- l1norm(beta)
  betaNormTraj <- oldNorm
  lassoCostTraj <- c()
  count <- 0
  lambda <- lambdaVec[j]
  
  while(l1norm(beta) >= 0 & count < 100){ # Lol i hard upperbounded #step by 100
    count <- count + 1
    oldNorm <- l1norm(beta)
    randInd <- floor(runif(1,1,6))
    randomTF <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
    randomTF[randInd] <- TRUE
    vecR <- as.matrix(r(Xtrain, Ytrain, beta))
    for (i in c(1:5)){ if (randomTF[i]) {
      print("hello") 
      print(i)
      betaCoordUpdate <- softThres(beta[i] + 1/N *Xtrain[,i]%*% vecR , lambda)
      beta[i] <- c(betaCoordUpdate)
    }
    }
    betaTraj <- rbind(betaTraj,beta)
    betaNormTraj <- append(betaNormTraj, l1norm(beta))
    lassoCostTraj <- append(lassoCostTraj, lassoCost(beta,Xtrain,Ytrain,lambda))
  }
  png(filename=paste("betaNorm,lambda =", lambda,".png"))
  plot(betaNormTraj, main = paste("betaNorm for lambda =", lambda),  xlab = "counter", ylab = "l1norm(beta)")
  dev.off()
  png(filename=paste("lassoCost, lambda =", lambda,".png"))
  plot(lassoCostTraj, main = paste("lassoCost for lambda =", lambda),  xlab = "counter", ylab = "lasso Cost")
  dev.off()
  lassoCostLambdaPath <- append(lassoCostLambdaPath,lassoCost(beta, Xtest, Ytest, lambda))
}
dev.off()
length(lassoCostLambdaPath)
plot(lassoCostLambdaPath) # find the cost on test sets along lambdas' path
lambdaCV <- lambdaVec[which.min(lassoCostLambdaPath)]



#######################################
### bootstrapping
####################################### 
for (k in c(1:10)){
  bootstrapTimes <- 30
  bootstrapBagSize <- 50
  lambda <- lambdaCV
  betaBootstrap <- c(0,0,0,0,0)
  for (i in c(1:bootstrapTimes)){
    bag <- floor(runif(bootstrapBagSize, 1, 51)) # probability that the bag has entires=51 is zero
    beta <- rep(0,5) # c(0,0,0,0,0)
    perm <- sample(c(1:50), 50, replace = TRUE, prob = NULL)
    Xtrain <- X[perm[1:40],]
    Ytrain <- Y[perm[1:40]]
    N <- dim(Xtrain)[1]
    Xtest <- X[perm[41:50],]
    Ytest <- Y[perm[41:50]]
    
    #betaTraj <- beta
    oldNorm <- l1norm(beta)
    betaNormTraj <- oldNorm
    lassoCostTraj <- c()
    count <- 0
    #lambda <- lambdaVec[j]
    
    while(l1norm(beta) >= 0 & count < 100){ # Lol i hard upperbounded #step by 100
      count <- count + 1
      oldNorm <- l1norm(beta)
      randInd <- floor(runif(1,1,6))
      randomTF <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
      randomTF[randInd] <- TRUE
      vecR <- as.matrix(r(Xtrain, Ytrain, beta))
      for (i in c(1:5)){ 
        if (randomTF[i]) {
          print("hello") 
          print(i)
          betaCoordUpdate <- softThres(beta[i] + 1/N *Xtrain[,i]%*% vecR , lambda)
          beta[i] <- c(betaCoordUpdate)
        }
      }
      #betaTraj <- rbind(betaTraj,beta)
      betaNormTraj <- append(betaNormTraj, l1norm(beta))
      lassoCostTraj <- append(lassoCostTraj, lassoCost(beta,Xtrain,Ytrain,lambda))
    }
    plot(betaNormTraj, main = paste("betaNorm for lambda =", lambda),  xlab = "counter", ylab = "l1norm(beta)")
    plot(lassoCostTraj, main = paste("lassoCost for lambda =", lambda),  xlab = "counter", ylab = "lasso Cost")
    betaBootstrap <- rbind(betaBootstrap,beta)
  }
  betaBootstrap <- betaBootstrap[-1,]
  betaStandardError <- apply(betaBootstrap, 2, sd) 
  # Reference http://stackoverflow.com/questions/18047896/column-standard-deviation-r
  # In apply(), the second argument, margin, "is a vector giving the subscripts ...
  # which the function will be applied over." 2 indicates columns 
  plot(betaStandardError, main = "beta Standard Error")
  write.matrix(betaTraj, file = paste("matrix", k, ".txt"), sep = " ")
  seBackup <- rbind(seBackup,c(30,lambdaMax, betaStandardError))
}
hist(seBackup[,3])

# CV how many times, lambdaMax
seBackup <- c(30, 0.34796, 0.02905589, 0.10109511, 0.06882610, 0.06479482, 0.02556103)
seBackup <- rbind(seBackup, c(30, 0.3995, 0.7555906, 1.9827836, 1.7068904, 1.5334808, 2.6828836))
seBackup <- rbind(seBackup, c(30, 0.3539, 5.653644,  9.475359, 11.753006, 3.375917, 4.350424))
seBackup <- rbind(seBackup, c(30, 0.3497, 0.08283931, 0.08443138, 0.20263155, 0.11515760, 0.14873422))
seBackup <- rbind(seBackup, c(30, 0.4651,  0.06320993, 0.05091544, 0.08786726, 0.01048634, 0.01892165))
