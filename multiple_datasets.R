library(MCMCpack)

source("logistic_regression.R")

logit_train_acc = c(0:10)
logit_test_acc = c(0:10)
mcmc_train_acc = c(0:10)
mcmc_test_acc = c(0:10)

for(i in 0:10)
{
  y = read.csv(paste("results/y_",i,".csv", sep=""), header=FALSE)
  y_test = read.csv(paste("results/y_test_",i,".csv", sep=""), header=FALSE)
  y = (y+1)/2
  y_test = (y_test+1)/2
  
  train = cbind(X,y)
  test = cbind(X_test,y_test)
  names(train) = c("x0","x1","x2","y")
  names(test) = c("x0","x1","x2","y")
  
  logReg = glm(y ~ x1 + x2, data = train, family = binomial)
  summary(logReg)
  w_logReg = coef(logReg)
  
  # Compute predicted probabilities on training data
  logPred = predict(logReg, type = "response")
  
  # Build a classification table to check accuracy on 
  # training set. 
  train.table = table(train$y, round(logPred))
  logit_train_acc[i+1] = sum(diag(train.table))/nrow(train)
  
  # We now do the same for the test set
  logPredTest = predict(logReg, newdata = test, type = "response")
  test.table <- table(test$y, round(logPredTest))
  test.table
  
  # Compute percentage correct (overall accuracy)
  logit_test_acc[i+1] = sum(diag(test.table))/nrow(test)
  
  # Run MCMC on dataset
  posterior = MCMClogit(y ~ x1 + x2, data = train, seed=0)
  w_mcmc = colMeans(posterior)[]
  V_mcmc = cov(posterior)
  colnames(V_mcmc) = c("x1","x2","x3")
  write.csv(as.matrix(w_mcmc), paste("results/w_mcmc_",i,".csv",sep=""), row.names=FALSE)
  write.csv(as.matrix(V_mcmc), paste("results/V_mcmc_",i,".csv",sep=""), row.names=FALSE)
}

df = data.frame(logit_train=logit_train_acc, logit_test=logit_test_acc)
write.csv(df, "results/accuracy_logit.csv")

