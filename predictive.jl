using DataFrames
using Gadfly
using Cairo
using Distributions
srand(1233)
include("mfvb_functions.jl")
include("vb_logit_fit.jl")
include("vb_pred.jl")
#include("vb_logit_fit_iter.jl")
M=10

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
y = convert(Array,readtable("data/example_coeff_y.csv", header=false))
y_test = convert(Array,readtable("data/example_coeff_y_test.csv", header=false))

accuracy_train=zeros(M+1)
accuracy_test=zeros(M+1)

for i = 0:M
	w=convert(Array,readtable("results/w_$i.csv", header=false))
	w=w[:,1]
	V=convert(Array,readtable("results/V_$i.csv", header=false))
	invV=convert(Array,readtable("results/invV_$i.csv", header=false))
	y=convert(Array,readtable("results/y_$i.csv", header=false))
	y=y[:,1]
	y_test=convert(Array,readtable("results/y_test_$i.csv", header=false))
	y_test=y_test[:,1]
	pred_train=vb_pred(X, w, V, invV)
	pred_test=vb_pred(X_test, w, V, invV)
	ypred_train=2*(pred_train .>0.5)-1
	ypred_test=2*(pred_test .>0.5)-1
	accuracy_train[i+1]=sum(ypred_train.==y)/length(y)
	accuracy_test[i+1]=sum(ypred_test.==y_test)/length(y_test)
end 
	
	accuracy=DataFrame(hcat(accuracy_train, accuracy_test))
	writetable("results/accuracy_MFVB.csv", accuracy, header=false)
	# writetable("results/y_test_$i.csv", y_test, header=false)
	# writetable("results/w_$i.csv", w, header=false)
	# writetable("results/V_$i.csv", V, header=false)
	# writetable("results/invV_$i.csv", invV, header=false)


