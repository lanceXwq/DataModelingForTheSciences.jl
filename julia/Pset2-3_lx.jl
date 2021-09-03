### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ dd274e79-5d97-4a60-8365-c0caefbaffba
begin
	import Pkg, Random
	Pkg.activate("")
	Random.seed!(1234)
	using Plots
end

# ╔═╡ d00c6ad6-af7c-11eb-3dbf-41fcda7e6cd4
md"
# Problem 3: Birth-death process
RNA is produced from DNA. 
RNA can also be degraded. 
Simulate this birth-death process for RNA (assuming a birth rate greater than the death rate). 
Plot the amount of RNA as a function of time until the value reaches steady state.
How is this steady state related to the birth and death rate?
"

# ╔═╡ c0420396-cf2f-4f07-b81a-6f3e1d05ca67
md"
## Solution:
"

# ╔═╡ 707f6c32-9311-47ae-9a71-67ddc1e816f7
md"
Start with specifying the birth rate and the death rate.
"

# ╔═╡ 2597b2a6-ff21-4d2c-a8bb-1fe09959d728
rate_birth = 100

# ╔═╡ d62088f7-06e9-4904-97ec-f6017a2b46ed
rate_death = 1

# ╔═╡ 6ed8fd8d-7ea9-4a8d-8441-ae57f27b70d2
N = 5000

# ╔═╡ 36eec2ba-f3d5-40f6-a9d0-26e6ce36c7a8
function gillespie(n_init::Int, rate_birth::Real, rate_death::Real, N::Int)
	n = zeros(Int64, N);
	t = zeros(Float64, N);
	n[1] = n_init;
	for idx = 2:N
		propensity = n[idx-1] * rate_death + rate_birth;
		t[idx] = t[idx-1] - log(rand())/propensity;
		if propensity * rand() > rate_birth
			n[idx] = n[idx-1] - 1;
		else
			n[idx] = n[idx-1] + 1;
		end
	end
	return n, t
end

# ╔═╡ 34b1bdf8-6147-425b-97a6-f4e31558f6be
(n, t) = gillespie(10, rate_birth, rate_death, 5000);

# ╔═╡ a666ea6c-5249-4e53-b036-50c7aa242f51
begin
	plot([0, t[end]], [rate_birth/rate_death, rate_birth/rate_death], linecolor = :red, label = "")
	plot!(t, n, linecolor = :steelblue, label = "")
end

# ╔═╡ Cell order:
# ╟─d00c6ad6-af7c-11eb-3dbf-41fcda7e6cd4
# ╟─c0420396-cf2f-4f07-b81a-6f3e1d05ca67
# ╠═dd274e79-5d97-4a60-8365-c0caefbaffba
# ╟─707f6c32-9311-47ae-9a71-67ddc1e816f7
# ╟─2597b2a6-ff21-4d2c-a8bb-1fe09959d728
# ╟─d62088f7-06e9-4904-97ec-f6017a2b46ed
# ╟─6ed8fd8d-7ea9-4a8d-8441-ae57f27b70d2
# ╠═36eec2ba-f3d5-40f6-a9d0-26e6ce36c7a8
# ╠═34b1bdf8-6147-425b-97a6-f4e31558f6be
# ╠═a666ea6c-5249-4e53-b036-50c7aa242f51
