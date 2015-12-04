function lam(xi)
	ans = [tanh(xi ./ 2) ./ (4 .* xi)]
	ans[isnan(ans)] = 0.125
	return ans
end