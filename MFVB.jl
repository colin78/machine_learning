using DataFrames
using Gadfly
using Cairo
using Distributions
srand(1233)
include("mfvb_functions.jl")
include("vb_logit_fit.jl")
#include("vb_logit_fit_iter.jl")
M=10

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))

#dic=Dict()
a0=rand(1,M)
b0=rand(1,M)
X_for_plot=0
X_for_plot2=0
for i=1:M
	#dic[i]=Dict()
	y=sampling(a0[i], b0[i], X)
	y_test=sampling(a0[i], b0[i], X_test)
	Matrix = vb_logit_fit(X,y)
	df = DataFrame(Matrix)
	w = DataFrame(w=df[:,1])
	V = df[:,2:4]
	invV=df[:, 5:7]
	writetable("results/y_$i.csv", y, header=false)
	writetable("results/y_test_$i.csv", y_test, header=false)
	writetable("results/w_$i.csv", w, header=false)
	writetable("results/V_$i.csv", V, header=false)
	writetable("results/invV_$i.csv", invV, header=false)
	#writetable("results/MFVBw.csv", w, header=false)
	#writetable("results/V-MFVB.csv", V, header=false)
	dic[i]["w"]=w
	dic[i]["V"]=V
	dic[i]["invV"]=V
	dic[i]["y"]=y
	dic[i]["y_test"]=y_test
	X_for_plot = DataFrame(x1=X[:,2],x2=X[:,3],y=y[:])
	X_for_plot2= DataFrame(x1=X_test[:,2],x2=X_test[:,3],y=y_test[:])
	
	equation(x1) = (-w[1] - w[2]*x1)/w[3]

	plot1 = plot(layer(X_for_plot, x="x1", y="x2", Geom.point, color = "y"),
				Scale.color_discrete_manual("red","blue"),
				layer(equation, -3, 3, Theme(default_color=color("orange"))), 
				Guide.xlabel("Stimulus"), Guide.ylabel("Response"), Guide.title("Dog Training"));
	draw(PDF("figures/train[$i].pdf", 6inch, 4inch), plot1)

	plot2 = plot(layer(X_for_plot2, x="x1", y="x2", Geom.point, color = "y"),
				Scale.color_discrete_manual("red","blue"),
				layer(equation, -3, 3, Theme(default_color=color("black"))), 
				Guide.xlabel("Stimulus"), Guide.ylabel("Response"), Guide.title("Dog Training"));
	draw(PDF("figures/test[$i].pdf", 6inch, 4inch), plot2)
end 


