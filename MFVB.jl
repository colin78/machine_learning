using DataFrames
using Gadfly
using Cairo
using Distributions
srand(1233)
include("mfvb_functions.jl")
include("vb_logit_fit.jl")
#include("vb_logit_fit_iter.jl")
M=100

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
#X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
#y = convert(Array,readtable("data/example_coeff_y.csv", header=false))
#y_test = convert(Array,readtable("data/example_coeff_y_test.csv", header=false))
dic=[]
a0=rand(1,M)
b0=rand(1,M)
for i=1:M
	dic[i]=[]
	y=sampling(a0[i], b0[i], X)
	Matrix = vb_logit_fit(X,y)
	df = DataFrame(Matrix)
	w = DataFrame(w=df[:,1])
	V = df[:,2:4]
	writetable("results/MFVBw.csv", w, header=false)
	writetable("results/V-MFVB.csv", V, header=false)

	files = ["w", "V", "invV"]
	fvars = Dict()
	for (n, f) in enumerate(files)
   		fvars["x_$(n)"] = f
	end
end 
y=sampling(a0, b0, X)
# X_for_plot = DataFrame(x1=X[:,2],x2=X[:,3],y=y[:])
# X_for_plot2 = DataFrame(x1=X_test[:,2],x2=X_test[:,3],y=y_test[:])

#Matrix_iter = vb_logit_fit_iter(X,y)
Matrix = vb_logit_fit(X,y)
df_iter = DataFrame(Matrix_iter)
df = DataFrame(Matrix)

w_iter = DataFrame(w_iter=df_iter[:,1])
w = DataFrame(w=df[:,1])
V_iter = df_iter[:,2:4]
V = df[:,2:4]

writetable("results/w_iter.csv", w_iter, header=false)
writetable("results/w.csv", w, header=false)
writetable("results/V_iter.csv", V_iter, header=false)
writetable("results/V.csv", V, header=false)