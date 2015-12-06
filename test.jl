using DataFrames

include("mfvb_functions.jl")

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
y = convert(Array,readtable("data/example_coeff_y.csv", header=false))
y_test = convert(Array,readtable("data/example_coeff_y_test.csv", header=false))
