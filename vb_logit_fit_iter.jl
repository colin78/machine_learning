using DataFrames

X = readtable("data/example_coeff_X.csv", header=false)
X_test = readtable("data/example_coeff_X_test.csv", header=false)
y = readtable("data/example_coeff_y.csv", header=false)
y_test = readtable("data/example_coeff_y_test.csv", header=false)

# function vb_logit_fit_iter(X,y)
N, D = size(X)
max_iter = 500

V = eye(D) / D
invV = eye(D) * D
logdetV = - D * log(D)
w = zeros(D, 1)

# return [X, y]
# end
