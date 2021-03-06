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
N,D=size(X)
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
Xt=vcat(X,X_test)
Nt,Dt=size(Xt)
dic=Dict()
a0=rand(1,M)
writetable("results/a0.csv", y, header=false)
b0=rand(1,M)
writetable("results/b0.csv", y, header=false)
X_for_plot=0
X_for_plot2=0
for i=1:M
	dic[i]=Dict()
	yt=sampling(a0[i], b0[i], Xt)
	y=yt[1:N,:]
	y_test=yt[N+1:Nt,:]
	Matrix = vb_logit_fit(X,y)
	df = DataFrame(Matrix)
	w=df[:,1]
	V = df[:,2:4]
	invV=df[:, 5:7]
	
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
				Guide.title("Separating Hyperplanes, Training Dataset$i"));
	draw(PDF("figures/train[$i].pdf", 6inch, 4inch), plot1)

	plot2 = plot(layer(X_for_plot2, x="x1", y="x2", Geom.point, color = "y"),
				Scale.color_discrete_manual("red","blue"),
				layer(equation, -3, 3, Theme(default_color=color("black"))), 
				Guide.title("Separating Hyperplanes, Testing Dataset$i"));
	draw(PDF("figures/test[$i].pdf", 6inch, 4inch), plot2)

	y=DataFrame(y)
	y_test=DataFrame(y_test)
	w=DataFrame(w=df[:,1])
	V = DataFrame(V)
	invV=DataFrame(invV)

	writetable("results/y_$i.csv", y, header=false)
	writetable("results/y_test_$i.csv", y_test, header=false)
	writetable("results/w_$i.csv", w, header=false)
	writetable("results/V_$i.csv", V, header=false)
	writetable("results/invV_$i.csv", invV, header=false)

end 


