using DataFrames

include("mfvb_functions.jl")

X = readtable("data/example_coeff_X.csv", header=false)
X_test = readtable("data/example_coeff_X_test.csv", header=false)
y = readtable("data/example_coeff_y.csv", header=false)
y_test = readtable("data/example_coeff_y_test.csv", header=false)

