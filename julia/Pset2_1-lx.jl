### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ d84527f2-cf85-407d-9d6c-e2e0c13c48b9
begin
	using LsqFit
	@. model_pdf(r, λ) = λ*exp(-r*λ);
	@. model_cdf(r, λ) = 1 - exp(-r*λ);
end

# ╔═╡ c1e65a92-f79a-4c11-8fae-b3b5ca31e8e5
md"
# Problem 1: Simulating stochastic processes
Imagine a process with an exponential waiting time between events as we have seen in class. 
We call this a birth process.\
a) Use the cdf of the exponential distribution to derive an expression for how to sample from an exponential distribution.\
b) Sample 100 exponential time intervals (you need to specify a rate, ``\lambda``, by hand). Use your generated data to construct the normalized cdf and pdf (you can choose arbitrary bins for the pdf histogram).\
c) Estimate the rate, ``\lambda``, from the mean of the exponential time of arrivals. 
Then estimate the mean by fitting the cdf and pdf. 
Which method for estimating ``\lambda`` seems more efficient?
"

# ╔═╡ 225ff9e3-bd16-4d72-8016-53c70ee6b095
md"
## Solution:
### a)
For an exponential distribution with the form $p(r)=\lambda e^{-\lambda r}$ for $r\geq 0$, the corresponding cdf, as shown in the example 1.13, is $C(r)=1-e^{-\lambda r}$. If follows that $C^{-1}(u)=-\frac{1}{\lambda}\log(1-u)$. According to the fundamental theorem of simulation, we first sample $u$ from $\text{Uniform}_{[0,1]}$, then $r$ is calculated using $r=C^{-1}(u)$.

"

# ╔═╡ d69b0313-ccc5-4407-9ae2-d68033adaf0d
function exprnd(λ, N)
    u = rand(N); # Generate u with some give dimension.
    return @. -1/λ * log(u) # Calculate r.
end;

# ╔═╡ 5534e093-7e22-4e6c-a813-681e2a687066
md"
Note that in the function ```exprnd``` ``u`` is used instead of ``(1-u)`` since they are equivalent.
"

# ╔═╡ c75b1010-7419-4e3a-be0f-77738dbb5417
md"
### b)
Now we draw 100 samples from the same exponential distribution and histogram the data to form the pdf.
"

# ╔═╡ 0ae23229-20be-4dd0-b75a-6a62c0e1dede
λ = 1

# ╔═╡ 6deb9836-71b1-4db4-8913-fbacb3863651
N = 100

# ╔═╡ fbdebc5d-406c-483c-8591-1ac39578ca78
r = exprnd(λ, N); # Call the functions we defined in a).

# ╔═╡ 047abd42-2095-4531-8779-221e18519555
begin
	using Plots, StatsBase, LinearAlgebra
	pdf = normalize(fit(Histogram, r), mode=:pdf)
	plot(pdf, label = "empirical")
	x = range(0, stop = maximum(pdf.edges...), length = 1000);
	plot!(x, λ*exp.(-λ .* x), label = "exact", linewidth = 2, linecolor = :red)
	xaxis!("r")
	yaxis!("pdf")
end

# ╔═╡ 4f0623c9-3a85-4110-bedb-56b58b124c21
md"Now the result can be tested via some visualization."

# ╔═╡ f8776b5f-43f5-411e-98a5-db596ab40013
md"
The cdf plots can either be made using the pdf histogram above, or they can be constructed directly. Here, we would show both ways.
"

# ╔═╡ 336f643e-aeb5-431b-be52-002a566e97e2
begin
	# Based on the pdf.
	width = pdf.edges[1][2] - pdf.edges[1][1];
	cdf_heights = cumsum(pdf.weights) * width;
	midpoints = (pdf.edges[1][1:end-1] + pdf.edges[1][2:end])/2;
	cdf = bar(midpoints, cdf_heights, label = "empirical 1", bar_width = width)
	y = range(0, stop = maximum(pdf.edges...), length = 1000);
	# Calculated directly.
	r_sorted = sort(r);
	plot!(r_sorted, (1:N)/N, label = "empirical 2", linewidth = 2, linecolor = :cyan)
	plot!(y, 1 .- exp.(-λ .* y), label = "exact", linewidth = 2, linecolor = :red)
	xaxis!("r")
	yaxis!("cdf")
end

# ╔═╡ 4690219e-b105-40f9-a574-fabf0f00925f
md"
### c)
"

# ╔═╡ 3f971df3-cbec-4a3c-bde2-8654e3528bc0
begin
	fit_pdf = curve_fit(model_pdf, midpoints, pdf.weights, [0.5]);
	fit_pdf.param
end

# ╔═╡ 152ad38e-bb79-4800-a1de-e9308714e306
begin
	fit_cdf1 = curve_fit(model_cdf, midpoints, cdf_heights, [0.5]);
	fit_cdf1.param # Get λ from fitting cdf 1.
end

# ╔═╡ d8d2e660-2d42-40ae-b095-8a50da826a9c
begin
	fit_cdf2 = curve_fit(model_cdf, r_sorted, (1:N)/N, [0.5]);
	fit_cdf2.param # Get λ from fitting cdf 2.
end

# ╔═╡ 47c7158f-186b-40ea-aba3-aaad0bcbc895
mean(r) # Get λ from calculating the mean of all samples.

# ╔═╡ Cell order:
# ╟─c1e65a92-f79a-4c11-8fae-b3b5ca31e8e5
# ╟─225ff9e3-bd16-4d72-8016-53c70ee6b095
# ╠═d69b0313-ccc5-4407-9ae2-d68033adaf0d
# ╠═5534e093-7e22-4e6c-a813-681e2a687066
# ╟─c75b1010-7419-4e3a-be0f-77738dbb5417
# ╟─0ae23229-20be-4dd0-b75a-6a62c0e1dede
# ╟─6deb9836-71b1-4db4-8913-fbacb3863651
# ╠═fbdebc5d-406c-483c-8591-1ac39578ca78
# ╟─4f0623c9-3a85-4110-bedb-56b58b124c21
# ╠═047abd42-2095-4531-8779-221e18519555
# ╟─f8776b5f-43f5-411e-98a5-db596ab40013
# ╠═336f643e-aeb5-431b-be52-002a566e97e2
# ╟─4690219e-b105-40f9-a574-fabf0f00925f
# ╠═d84527f2-cf85-407d-9d6c-e2e0c13c48b9
# ╠═3f971df3-cbec-4a3c-bde2-8654e3528bc0
# ╠═152ad38e-bb79-4800-a1de-e9308714e306
# ╠═d8d2e660-2d42-40ae-b095-8a50da826a9c
# ╠═47c7158f-186b-40ea-aba3-aaad0bcbc895
