include("mfvb_functions.jl")

function vb_logit_fit(X,y)
    # hyperprior parameters
    a0 = 1e-2
    b0 = 1e-4
    # pre-compute some constants
    N, D = size(X)
    max_iter = 500
    an = a0 + 0.5 * D;    gammaln_an_an = lgamma(an) + an
    t = 0.5 * sum( X.*y,1)'

    # start first iteration kind of here, with xi = 0 -> lam_xi = 1/8
    lam_xi = ones(N, 1) / 8
    E_a = a0 / b0
    invV = E_a * eye(D) + 2 * X' *  (X.*lam_xi)
    V = inv(invV)
    w = V * t
    bn = b0 + 0.5 * (w' * w + trace(V))
    bn=bn[1]
    L_last = - N * log(2) + 0.5 * (w' * invV * w - logdet(invV)) - an / bn * b0 - an * log(bn) + gammaln_an_an
    L_last = L_last[1]
    L = 0

    #update xi, bn, (V, w) iteratively
    for i = 1:max_iter
        #update xi by EM-algorithm
        xi = sqrt(sum(X .* (X * (V + w * w')), 2))
        xi = convert(Array{Float64,1},xi[:,1])
        lam_xi = lam(xi)

        #update posterior parameters of a based on xi
        bn = b0 + 0.5 * (w' * w + trace(V))
        E_a = an / bn[1]

        #recompute posterior parameters of w
        invV = E_a * eye(D) + 2 * X' *  (X.*lam_xi)
        V = inv(invV)
        logdetV = - logdet(invV)
        w = V * t

        #variational bound, ignoring constant terms for now
        L = - sum(log(1 + exp(- xi))) + sum(lam_xi .* xi .^ 2) + 0.5 * (w' * invV * w + logdetV - sum(xi))- E_a * b0 
        - an * log(bn) + gammaln_an_an
        L=L[1]

         # either stop if variational bound grows or change is < 0.001%
         # HACK ALARM: theoretically, the bound should never grow, and it doing
         # so points to numerical instabilities. As it seems, these start to
         # occur close to the optimal bound, which already points to a good
         # approximation.
        if (L_last > L) | (abs(L_last - L) < abs(0.00001 * L))
            break
        end
        L_last = L

        if i == max_iter
            println("Bayesian logistic regression reached maximum number of iterations")
        end
    end

    #add constant terms to variational bound
    L = L - lgamma(a0) + a0 * log(b0)

    return [w V]
end


