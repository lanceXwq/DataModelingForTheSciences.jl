### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ a3ac2089-cea2-40b8-bf2e-edc212c611e7
using LaTeXStrings, Plots, StatsBase, LinearAlgebra

# ╔═╡ 206cb0d2-a843-11eb-0b72-f942682593b2
md"
# Problem 1: EM algorithm
Implement the EM algorithm for the simple case of a mixture of two Gaussians in 1D with identical standard deviation. 
In other words, begin by generating synthetic data according to this model: ``\dfrac{\pi_1}{\sqrt{2\pi\sigma^2}} \exp\left[-\dfrac{(x-\mu_1)^2}{2\sigma^2}\right]+\dfrac{\pi_2}{\sqrt{2\pi\sigma^2}} \exp\left[-\dfrac{(x-\mu_2)^2}{2\sigma^2}\right]``, 
where you have pre-specified by hand the parameters ``\pi_1,\mu_1,\mu_2,\sigma``
Then, implement an EM in order to learn the parameters ``\pi_1,\mu_1,\mu_2,\sigma``. 
Compare your parameter estimates from EM to the theoretical values you used to generate your data.
"

# ╔═╡ 70b8d1da-bac4-4d65-912a-bce179e632ca
md"
## Solution:
"

# ╔═╡ 32ef202f-9854-492f-a5cc-67de57f7dd9a
md"
### Synthetic data
Let's start from specifiying the values of ``\pi_1,\mu_1,\mu_2,\sigma``, as well as the number of data points, ``N``.
"

# ╔═╡ cefb0e0e-4155-402d-a2cb-9a7e38105d4f
π₁ = 0.6

# ╔═╡ 31e09874-f041-408b-b2f7-a6e25a0aa45d
μ=[3, 5]

# ╔═╡ a86b5d8a-99d6-4f9b-9ba8-36992e14ffea
σ = 1

# ╔═╡ 034b2209-3618-4423-9097-3af3d7299055
N = 1000

# ╔═╡ dd69e9dc-9b6b-42eb-86b5-e1f751a711d0
typeof(μ)

# ╔═╡ 1e89f4ed-be1a-4408-ac33-0ca58656c904
gaussian(y, μ::Real, σ::Real) = @. 1 / √(2*π*σ^2) * exp(-(y-μ)^2/(2*σ^2))

# ╔═╡ 1307292d-c51c-4ce8-abd0-7de9d9161a00
GMM(y, pi::Real, μ, σ::Real) = @. pi * gaussian(y, μ[1], σ) + (1 - pi) * gaussian(y, μ[2], σ)

# ╔═╡ d785ad01-7d14-43f1-9685-12ace5f631d7
md"
Generating synthetic data, in this case, can be divided into two steps
1. sample from a categorical distribution for the state,
2. sample from a Gaussian distribution for the observtion.
Since we only have two possible states in this problem, we can write a simple and fast categorical distribution sampler called ```randc```.
"

# ╔═╡ 2595ced6-2285-4a4b-ab6a-9e0f1ad193b1
randc(p::Real, N::Integer) = (rand(N) .> p) .+ 1

# ╔═╡ b4039295-3ce1-42c1-a43d-eb8be79d0105
md"
Then the state of each data point is simply ```randc(π₁, N)```.
"

# ╔═╡ 9edea97a-bbfb-454d-9430-c2d4a2609d28
s = randc(π₁, N);

# ╔═╡ 5dffde03-8917-4a33-a939-0b2889c45af2
md"
The observations are obtaind via adding noise sampled from the Gaussiam distributions.
"

# ╔═╡ b8fb8905-a288-459e-828f-6765f7883459
y = μ[s] .+ σ .* randn(N);

# ╔═╡ c323acb0-8b89-49bc-a93a-efc641e76909
begin
	pdf = normalize(fit(Histogram, y), mode=:pdf)
	plot(pdf, label = "histogram")
	x = range(0, stop = maximum(pdf.edges...), length = 1000);
	plot!(x, GMM(x, π₁, μ, σ), label = "exact", linewidth = 2, linecolor = :red)
	xaxis!("data")
	yaxis!("pdf")
end

# ╔═╡ 0badca44-d32b-4771-b6f8-e66624551add
md"
## EM algorithm
Most of the derivation needed has already been done in the textbook, here, recall that we have defined
```math
\gamma_{1n}^{old}=\frac{\pi_1^{old}\exp\left[\frac{-\left(y_n-\mu_1\right)^{2}}{2 \sigma^{2}}\right]}{\pi_1^{old}\exp\left[\frac{-\left(y_n-\mu_1\right)^{2}}{2 \sigma^{2}}\right]+\pi_2^{old}\exp\left[\frac{-\left(y_n-\mu_2\right)^{2}}{2 \sigma^{2}}\right]},
```
and
```math
\gamma_{2n}^{old}=1-\gamma_{1n}^{old}.
```
"

# ╔═╡ d6b05d9e-a1c4-431a-b9c1-a65fa961e6f6
function gamma(y, pi::Real, μ::Vector{<:Real}, σ::Real) 
	return @. 1 / (1 + (1 - pi) / pi * exp(((y-μ[1])^2-(y-μ[2])^2)/(2*σ^2)))
end

# ╔═╡ 5fb43039-584b-4f11-957f-84edc3ae996c
md"
Then the first four parameters are updated as
```math
\pi_1^{new}=\frac{\sum_{n=1}^N \gamma_{1n}^{old}}{\sum_{n=1}^N \left(\gamma_{1n}^{old}+\gamma_{2n}^{old}\right)}=\frac{1}{N}\sum_{n=1}^N \gamma_{1n}^{old},~~~~\pi_2^{new}=1-\pi_1^{new},
```
and
```math
\mu_1^{new}=\frac{\sum_{n=1}^N \gamma_{1n}^{old}y_n}{\sum_{n=1}^N \gamma_{1n}^{old}},
~~~~\mu_2^{new}=\frac{\sum_{n=1}^N \gamma_{2n}^{old}y_n}{\sum_{n=1}^N \gamma_{2n}^{old}}.
```
Learning ``\sigma`` requires a bit extra derivation but it is simple, and the result is
```math
\left(\sigma^2\right)^{new}=\frac{1}{N}\sum_{n=1}^N \left[\gamma_{1n}^{old}\left(y_n-\mu_1^{new}\right)^2+\gamma_{2n}^{old}\left(y_n-\mu_2^{new}\right)^2\right].
```
"

# ╔═╡ 7dfd8a5b-6b00-4a41-a4de-a5b927a3e4dd
md"
Now we start coding by specifying the number of iterations and guessing some initial values for the parameters.
"

# ╔═╡ e43a0915-a2a7-4cfa-86ee-867b29ba08fb
NIter = 1000

# ╔═╡ f5814ae5-59f1-46a6-be7e-c49aa6c5983c
π₁_old = zeros(Float64, NIter, 1);

# ╔═╡ 8d8321a0-0126-4445-a825-b76025ab3832
π₁_old[1] = 0.3;

# ╔═╡ 35fb2f7f-06f4-48c5-a6cc-a87bd54eca0a
μ_old = zeros(Float64, NIter, 2);

# ╔═╡ 12d45a86-776f-45dc-b741-86ff88f6295e
μ_old[1, :] = [2, 8];

# ╔═╡ 94c587ea-74b2-402e-b72c-aa24c55a7287
σ_old = zeros(Float64, NIter, 1);

# ╔═╡ 5a5ed5c9-b82b-41fd-b810-ef1df4f0ad8d
σ_old[1] = 2;

# ╔═╡ ad850209-4383-4754-8137-e083a6753549
for idx = 1:NIter - 1
	gamma1_old = gamma(y, π₁_old[idx], μ_old[idx, :], σ_old[idx])
	gamma2_old = 1 .- gamma1_old
	gamma1_sum = sum(gamma1_old)
	gamma1_y_sum = sum(gamma1_old .* y)
	gamma2_y_sum = sum(gamma2_old .* y)
	π₁_old[idx+1] = gamma1_sum / N
	μ_old[idx+1, 1] = gamma1_y_sum / gamma1_sum
	μ_old[idx+1, 2] = gamma2_y_sum / (N - gamma1_sum)
	σ_old[idx+1] = sum(gamma1_old .* (y .- μ_old[idx+1, 1]).^2 + gamma2_old .* (y .- μ_old[idx+1, 2]).^2)
	σ_old[idx+1] = sqrt(σ_old[idx+1] / N)
end

# ╔═╡ c60dd468-234b-4dec-a1bb-c0e87d4f4b14
md"
## Result visualization
"

# ╔═╡ 1f2c55fc-1414-4274-b2f3-4c346f9eef39
begin
	plot([1, NIter], [π₁, π₁], label = "ground truth", linecolor = :red, line = (:dash, 4))
	plot!(π₁_old, label = "EM", linecolor = :steelblue)
	xaxis!("Iteration")
	yaxis!(L"\pi_1")
end

# ╔═╡ 932aaef6-5f49-4eea-aab9-90259998e72d
begin
	plot([1, NIter], [μ[1], μ[1]], label = "ground truth", linecolor = :red, line = (:dash, 4))
	plot!([1, NIter], [μ[2], μ[2]], label = "", linecolor = :red, line = (:dash, 4))
	plot!(μ_old[:, 1], label = "EM", linecolor = :steelblue)
	plot!(μ_old[:, 2], label = "", linecolor = :steelblue)
	xaxis!("Iteration")
	yaxis!(L"\mu")
end

# ╔═╡ c2520042-8977-48e7-9811-0fd435958053
begin
	plot([1, NIter], [σ, σ], label = "ground truth", linecolor = :red, line = (:dash, 4))
	plot!(σ_old, label = "EM", linecolor = :steelblue)
	xaxis!("Iteration")
	yaxis!(L"\pi_1")
end

# ╔═╡ Cell order:
# ╟─206cb0d2-a843-11eb-0b72-f942682593b2
# ╟─70b8d1da-bac4-4d65-912a-bce179e632ca
# ╠═a3ac2089-cea2-40b8-bf2e-edc212c611e7
# ╟─32ef202f-9854-492f-a5cc-67de57f7dd9a
# ╟─cefb0e0e-4155-402d-a2cb-9a7e38105d4f
# ╟─31e09874-f041-408b-b2f7-a6e25a0aa45d
# ╟─a86b5d8a-99d6-4f9b-9ba8-36992e14ffea
# ╟─034b2209-3618-4423-9097-3af3d7299055
# ╠═dd69e9dc-9b6b-42eb-86b5-e1f751a711d0
# ╠═1e89f4ed-be1a-4408-ac33-0ca58656c904
# ╠═1307292d-c51c-4ce8-abd0-7de9d9161a00
# ╟─d785ad01-7d14-43f1-9685-12ace5f631d7
# ╠═2595ced6-2285-4a4b-ab6a-9e0f1ad193b1
# ╟─b4039295-3ce1-42c1-a43d-eb8be79d0105
# ╠═9edea97a-bbfb-454d-9430-c2d4a2609d28
# ╟─5dffde03-8917-4a33-a939-0b2889c45af2
# ╠═b8fb8905-a288-459e-828f-6765f7883459
# ╠═c323acb0-8b89-49bc-a93a-efc641e76909
# ╟─0badca44-d32b-4771-b6f8-e66624551add
# ╠═d6b05d9e-a1c4-431a-b9c1-a65fa961e6f6
# ╟─5fb43039-584b-4f11-957f-84edc3ae996c
# ╟─7dfd8a5b-6b00-4a41-a4de-a5b927a3e4dd
# ╟─e43a0915-a2a7-4cfa-86ee-867b29ba08fb
# ╠═f5814ae5-59f1-46a6-be7e-c49aa6c5983c
# ╠═8d8321a0-0126-4445-a825-b76025ab3832
# ╠═35fb2f7f-06f4-48c5-a6cc-a87bd54eca0a
# ╠═12d45a86-776f-45dc-b741-86ff88f6295e
# ╠═94c587ea-74b2-402e-b72c-aa24c55a7287
# ╠═5a5ed5c9-b82b-41fd-b810-ef1df4f0ad8d
# ╠═ad850209-4383-4754-8137-e083a6753549
# ╟─c60dd468-234b-4dec-a1bb-c0e87d4f4b14
# ╠═1f2c55fc-1414-4274-b2f3-4c346f9eef39
# ╠═932aaef6-5f49-4eea-aab9-90259998e72d
# ╠═c2520042-8977-48e7-9811-0fd435958053
