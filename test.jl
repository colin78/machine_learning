using DataFrames
using Gadfly
using Cairo
using Distributions
srand(1233)
include("mfvb_functions.jl")
include("vb_logit_fit.jl")
include("vb_logit_fit_iter.jl")

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
y = convert(Array,readtable("data/example_coeff_y.csv", header=false))
y_test = convert(Array,readtable("data/example_coeff_y_test.csv", header=false))

X_for_plot = DataFrame(x1=X[:,2],x2=X[:,3],y=y[:])
X_for_plot2 = DataFrame(x1=X_test[:,2],x2=X_test[:,3],y=y_test[:])

Matrix_iter = vb_logit_fit_iter(X,y)
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


equation(x1) = (-w_iter[1] - w_iter[2]*x1)/w_iter[3]
equation2(x1) = (-w[1] - w[2]*x1)/w[3]

# plot1 = plot(layer(X_for_plot, x="x1", y="x2", Geom.point, color = "y"),
# 			Scale.color_discrete_manual("red","blue"),
# 			layer(equation, -3, 3, Theme(default_color=color("green"))),
# 			layer(equation2, -3, 3, Theme(default_color=color("orange"))));
# draw(PDF("figures/train.pdf", 6inch, 4inch), plot1)

# plot2 = plot(layer(X_for_plot2, x="x1", y="x2", Geom.point, color = "y"),
# 			Scale.color_discrete_manual("red","blue"),
# 			layer(equation2, -3, 3, Theme(default_color=color("black"))));
# draw(PDF("figures/test.pdf", 6inch, 4inch), plot2)