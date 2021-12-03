# 2(b)
library(glmnet)
library(pracma)
library(matlib)
# simulate a dataset with 100 samples, 200 features (random draws from U(-1,1)) using equation (1) with epsilon~N(0,1). 
# Obtain beta directly from equation (3). 
set.seed(123)

#define matrix of predictor variables
x <- replicate(200, runif(n=100,min=-1,max=1))
beta=rep(1,length=200)


#define response variable
epsilon<- rnorm(n=100,mean=0,sd=1)
y=x%*%beta+epsilon


# need to find the optimal lambda

#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 0)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda

start.time <- Sys.time()
beta_hat = (inv(t(x)%*%x+best_lambda*I))%*%(t(x)%*%y) 
beta_hat
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# 2(c)use adaptive gradient descent  (gradDescent package in R) to obtain  ?? 

library(gradDescent)
df=data.frame(x,y)
splitdf=splitData(df,dataTrainRate = 0.75,seed=111)

start.time <- Sys.time()
ADAGRADmodel <-ADAGRAD(splitdf$dataTrain) 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#show result 
#print(ADAGRADmodel)
t(ADAGRADmodel)
