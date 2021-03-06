include("mfvb_functions.jl")

function vb_logit_fit_iter(X,y)
	N, D = size(X)
	max_iter = 500

	V = eye(D) / D
	invV = eye(D) * D
	logdetV = - D * log(D)
	w = zeros(D, 1)

	for n = 1:N
	    xn = X[n,:]'

	    # precompute values
	    Vx = V * xn
	    VxVx = Vx * Vx'
	    c = (xn' * Vx)[1]
	    xx = xn * xn'
	    t_w = invV * w + 0.5 * y[n] * xn

	    # start iteration at xi = 0, lam_xi = 1/8
	    V_xi = V - VxVx / (4 + c)
	    invV_xi = invV + xx / 4
	    logdetV_xi = logdetV - log(1 + c / 4)
	    w = V_xi * t_w
	    
	    L_last = 0.5 * (logdetV_xi + w' * invV_xi * w)[1] - log(2)

	    for i = 1:max_iter
	        # update xi by EM algorithm
	        xi = sqrt((xn' * (V_xi + w * w') * xn)[1])
	        lam_xi = lam(xi)[1]

	        # Sherman-Morrison formula and Matrix determinant lemma
	        V_xi = V - (2 * lam_xi / (1 + 2 * lam_xi * c))[1] * VxVx
	        invV_xi = invV + 2 * lam_xi * xx
	        logdetV_xi = logdetV - log(1 + 2 * lam_xi * c)
	        w = V_xi * t_w

	        L = (0.5 * (logdetV_xi + w' * invV_xi * w - xi) - log(1 + exp(- xi)) + lam_xi * xi^2)[1]

	        # variational bound must grow!
	        if L_last > L
	            println("Last bound %6.6f, current bound %6.6f\n", L_last, L)
	            println("Variational bound should not reduce")
	        end
	        # stop if change in variation bound is < 0.001%
	        if abs(L_last - L) < abs(0.00001 * L)
	            break
	        end
	        L_last = L
	        if i == max_iter
		        println("Bayesian logistic regression reached maximum number of iterations.")
		    end
	    end
	    
	    V = V_xi
	    invV = invV_xi
	    logdetV = logdetV_xi
	end

	return [w V invV]
end

