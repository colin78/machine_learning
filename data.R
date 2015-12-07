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