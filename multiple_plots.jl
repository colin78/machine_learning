using DataFrames
using Gadfly
using Cairo
using Distributions

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))

for i=0:10
	w_mfvb = readtable("results/w_$i.csv", header=false)[:,1]
	w_mcmc = readtable("results/w_mcmc_$i.csv", header=true)[:,1]
	w_logit = readtable("results/w_logit_$i.csv", header=true)[:,1]
	y = readtable("results/y_$i.csv", header=false)[:,1]
	y_test = readtable("results/y_test_$i.csv", header=false)[:,1]

	equation_mfvb(x1) = (-w_mfvb[1] - w_mfvb[2]*x1)/w_mfvb[3]
	equation_mcmc(x1) = (-w_mcmc[1] - w_mcmc[2]*x1)/w_mcmc[3]
	equation_logit(x1) = (-w_logit[1] - w_logit[2]*x1)/w_logit[3]

	X_for_plot = DataFrame(x1=X[:,2],x2=X[:,3],y=y[:])
	X_for_plot_test = DataFrame(x1=X_test[:,2],x2=X_test[:,3],y=y_test[:])

	plot1 = plot(layer(X_for_plot, x="x1", y="x2", Geom.point, color = "y"),
				Scale.color_discrete_manual("red","blue"),
				layer(equation_mfvb, -3, 3, Theme(default_color=color("orange"))),
				layer(equation_mcmc, -3, 3, Theme(default_color=color("green"))),
				layer(equation_logit, -3, 3, Theme(default_color=color("black"))));
	draw(PDF("figures/train_final_$i.pdf", 6inch, 4inch), plot1)

	plot_test = plot(layer(X_for_plot_test, x="x1", y="x2", Geom.point, color = "y_test"),
				Scale.color_discrete_manual("red","blue"),
				layer(equation_mfvb, -3, 3, Theme(default_color=color("orange"))),
				layer(equation_mcmc, -3, 3, Theme(default_color=color("green"))),
				layer(equation_logit, -3, 3, Theme(default_color=color("black"))));
	draw(PDF("figures/test_final_$i.pdf", 6inch, 4inch), plot_test)
end