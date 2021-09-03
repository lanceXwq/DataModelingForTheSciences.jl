### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 2929689f-a9fd-4210-91e7-26505c9b8c9a
begin
	import Pkg, Random
	Pkg.activate("")
	Random.seed!(1234)
	using SpecialFunctions, Statistics, Plots, StatsBase, LinearAlgebra
end

# ╔═╡ 98a4942c-fca5-4d42-b043-df24e76aea23
md"
# Problem 2: Poisson processes
a) Devise a sampling scheme to sample from the Poisson distribution and code it up (i.e. don’t use the built in packages).\
b) Generate 1000 Poisson random variables from the Poisson distribution. 
Compute the ratio of mean and variance after every random variable is collected (i.e., using only the first and second data points; then using the first, second and third, etc...).
Plot this ratio from 2 to 1000 (where, just to be clear, for 1000 you would be using all data points to compute the mean and variance). To what value does this ratio converge?
"

# ╔═╡ 2fcf069a-729a-4918-be1c-9e98e8060bdd
md"
## Solution:
### a)
Naively, sampling from the Poisson distribution means we need to solve for the smallest non-negative integer ``N`` such that 
```math
\sum_{n=0}^N \frac{e^{-\lambda}\lambda^n}{n!}\geq u,
```
where ``u`` is sampled from the uniform distribution between 0 and 1.
However, this inequality is not trivial to solve, so we may want to have a workaround. 
Consider that the Poisson distribution expresses the probability of a given number of events occurring in a fixed interval of time or space if these events occur with a known constant mean rate and independently of the time since the last event. 
In other words, we can keep sampling individual events from an identical exponential distribution until it hits the length of the interval.
Then the problem becomes finding the largest ``N`` such that
```math
\sum_{n=1}^N -\frac{1}{\lambda}\ln{u_n}=-\frac{1}{\lambda}\ln\left(\prod_{n=1}^N u_n\right)\leq 1,
```
or equivalently,
```math
\prod_{n=1}^N u_n\leq e^{-\lambda}.
```
"

# ╔═╡ 90905632-a925-4d08-93d5-7ad8d6dc54d1
function poissrnd(λ::Real, N::Int)
	n = zeros(Int, N)
	lim = exp(-λ)
	for idx = 1:N
		u = rand()
		while u > lim
			u = u * rand()
			n[idx] += 1
		end
	end
    return n
end

# ╔═╡ fb6513cb-1f6d-4bbe-95de-45e47b0d2d44
md"
### b)
"

# ╔═╡ a0352d17-18e9-43c4-8f75-0f3e4ce40358
r = poissrnd(20, 1000);

# ╔═╡ 0151c22f-6ce7-4962-8978-7e2f5478824d
function poisspdf(n::Int, λ::Real)
	return exp(n * log(λ) - λ - logfactorial(n))
end

# ╔═╡ 5b55e3af-9cb1-493c-a57a-5a7a9b4383ee
begin
	pdf = normalize(fit(Histogram, r), mode=:pdf)
	plot(pdf, label = "empirical")
	x = collect(0:Int(round(maximum(pdf.edges...))))
	plot!(x, poisspdf.(x,20), label = "exact", linewidth = 2, linecolor = :red)
	xaxis!("r")
	yaxis!("pdf")
end

# ╔═╡ b61878e3-fe11-470e-8cca-62d35a39d7e9
begin
	mean_val = zeros(Float64, 999);
	var_val = zeros(Float64, 999);
	for idx = 2:1000
		mean_val[idx-1] = mean(r[1:idx]);
		var_val[idx-1] = var(r[1:idx]);
	end
	plot([1,1000], [1,1], label = "", linecolor = :red)
	plot!(2:1000, mean_val ./ var_val, label = "", linecolor = :steelblue)
end

# ╔═╡ Cell order:
# ╟─98a4942c-fca5-4d42-b043-df24e76aea23
# ╠═2929689f-a9fd-4210-91e7-26505c9b8c9a
# ╟─2fcf069a-729a-4918-be1c-9e98e8060bdd
# ╠═90905632-a925-4d08-93d5-7ad8d6dc54d1
# ╟─fb6513cb-1f6d-4bbe-95de-45e47b0d2d44
# ╠═a0352d17-18e9-43c4-8f75-0f3e4ce40358
# ╠═0151c22f-6ce7-4962-8978-7e2f5478824d
# ╠═5b55e3af-9cb1-493c-a57a-5a7a9b4383ee
# ╠═b61878e3-fe11-470e-8cca-62d35a39d7e9
