### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 7a51cb04-9e41-11ec-2c7c-39e816fd07c3
begin
    import Pkg
    Pkg.activate(Base.current_project())
    using GLMakie
    using Random
    using Distributions
    Random.seed!(1234)
end

# ╔═╡ 7a6605b9-07bf-445c-b6ee-f916304c215f
target = truncated(Normal(3, 2), lower=0, upper=5)

# ╔═╡ 390eb47e-45b2-4d0b-a4ee-125400a3fb3d
function get_accep_ratio(
    target::UnivariateDistribution,
    Q_old::UnivariateDistribution,
    Q_prop::UnivariateDistribution,
    old::Real,
    prop::Real,
)
    return min(
        1,
        pdf(target, prop) * pdf(Q_prop, old) / (pdf(target, old) * pdf(Q_old, prop)),
    )
end

# ╔═╡ e8c264eb-c4ec-43a0-954d-92ef7ea03b61
r = zeros(1000)

# ╔═╡ 40b20e7d-dc75-4186-80d9-bc28d46eab97
r[1] = 1

# ╔═╡ 11040e0a-85f4-40b4-95f8-7f2424f16467
for idx = 2:1000
    Q_old = Normal(r[idx-1], 2)
    prop = rand(Q_old)
    Q_prop = Normal(prop, 2)
    if get_accep_ratio(target, Q_old, Q_prop, r[idx-1], prop) > rand()
        r[idx] = prop
    else
        r[idx] = r[idx-1]
    end
end

# ╔═╡ d2c60d67-b5ca-42a0-bf76-41c5e294c3e9
begin
    x = -0:0.1:5
    hist(r, normalization=:pdf)
    lines!(x, pdf.(target, x), color=:red)
    current_figure()
end

# ╔═╡ Cell order:
# ╠═7a51cb04-9e41-11ec-2c7c-39e816fd07c3
# ╠═7a6605b9-07bf-445c-b6ee-f916304c215f
# ╠═390eb47e-45b2-4d0b-a4ee-125400a3fb3d
# ╠═e8c264eb-c4ec-43a0-954d-92ef7ea03b61
# ╠═40b20e7d-dc75-4186-80d9-bc28d46eab97
# ╠═11040e0a-85f4-40b4-95f8-7f2424f16467
# ╠═d2c60d67-b5ca-42a0-bf76-41c5e294c3e9
