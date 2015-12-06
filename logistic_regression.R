X = read.csv("data/example_coeff_X.csv", header=FALSE)
y = read.csv("data/example_coeff_y.csv", header=FALSE)
X_test = read.csv("data/example_coeff_X_test.csv", header=FALSE)
y_test = read.csv("data/example_coeff_y_test.csv", header=FALSE)

y = (y+1)/2
y_test = (y_test+1)/2

train = cbind(X,y)
test = cbind(X_test,y_test)

names(train) = c("x0","x1","x2","y")
names(test) = c("x0","x1","x2","y")

logReg = glm(y ~ x1 + x2, data = train, family = binomial)
summary(logReg)

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