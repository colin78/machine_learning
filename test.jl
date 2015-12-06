using DataFrames
using Gadfly
using Cairo

include("mfvb_functions.jl")

X = convert(Array,readtable("data/example_coeff_X.csv", header=false))
X_test = convert(Array,readtable("data/example_coeff_X_test.csv", header=false))
y = convert(Array,readtable("data/example_coeff_y.csv", header=false))
y_test = convert(Array,readtable("data/example_coeff_y_test.csv", header=false))

X_for_plot = DataFrame(x1=X[:,2],x2=X[:,3],y=y[:])
X_for_plot2 = DataFrame(x1=X_test[:,2],x2=X_test[:,3],y=y_test[:])


# for i = 1:n
# 	X_for_plot[i,:mip_cl] = indmax(mip_matrix[i,:])
# 	X_for_plot[i,:k_means_cl] = indmax(k_means_matrix[i,:])
# end

plot1 = plot(layer(X_for_plot, x="x1", y="x2", Geom.point, color = "y"),
			Scale.color_discrete_manual("red","blue"));
draw(PDF("figures/train.pdf", 6inch, 4inch), plot1)

plot1 = plot(layer(X_for_plot, x="x1", y="x2", Geom.point, color = "y"),
			Scale.color_discrete_manual("red","blue"));
draw(PDF("figures/train.pdf", 6inch, 4inch), plot1)