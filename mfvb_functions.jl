function lam(xi)
	ans = [tanh(xi ./ 2) ./ (4 .* xi)]
	ans[isnan(ans)] = 0.125
	return ans
end

function sampling(a0, b0,X)
	N,D=size(X)
	d1=Gamma(a0, b0)
	alpha=rand(d1)
	d2=MvNormal(eye(N)./alpha)
	w=rand(d2)
	y=zeros(N,1)
	for i=1:N
		x=X[i,:]
		p=1./(1+exp(-x*w))
		a=rand()
		if a < p[1] 
			y[i]=1
		else 
			y[i]=-1
		end 
	end 
	return [y]
end 