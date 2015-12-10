source("data.R")

logReg = glm(y ~ x1 + x2, data = train, family = binomial)
summary(logReg)
w_logReg = coef(logReg)

# Compute predicted probabilities on training data
logPred = predict(logReg, type = "response")

# Build a classification table to check accuracy on 
# training set. 
table(train$y, round(logPred))

# We now do the same for the test set
logPredTest = predict(logReg, newdata = test, type = "response")
test.table <- table(test$y, round(logPredTest))
test.table

# Compute percentage correct (overall accuracy)
sum(diag(test.table))/nrow(test)
