using GLMakie
using Random
using ColorSchemes

Random.seed!(1)

function get_brownian_traj(init_x, D, t, num_of_particle)
    Δt = diff(t)
    pushfirst!(Δt, 0)
    return init_x .+ cumsum(.√(2 * D * Δt) .* randn(length(Δt), num_of_particle), dims=1)
end

normpdf(x, μ, σ) = 1 / √(2π * σ^2) * exp.(-(x .- μ) .^ 2 ./ (2 * σ^2))

start_time = 0
end_time = 1
Δt = 0.01
t = start_time:Δt:end_time

init_x = 0

num_of_particle = 500

D = [0, 1, 4]

fig = Figure()

axes = [Axis(fig[1, 1]), Axis(fig[2, 1]), Axis(fig[3, 1]), Axis(fig[1, 2]), Axis(fig[2, 2]), Axis(fig[3, 2])]
linkaxes!(axes[1], axes[2], axes[3])
linkxaxes!(axes[4], axes[5], axes[6])

x1 = get_brownian_traj(init_x, D[1], t, num_of_particle)
for x in eachcol(x1)
    lines!(axes[1], t, x)
end

x2 = get_brownian_traj(init_x, D[2], t, num_of_particle)
for x in eachcol(x2)
    lines!(axes[2], t, x)
end

x3 = get_brownian_traj(init_x, D[3], t, num_of_particle)
for x in eachcol(x3)
    lines!(axes[3], t, x)
end

xlims!(axes[3], t[1], t[end])

slider_t = Slider(fig[4, 1], range=1:length(t), startvalue=0)

slice_pos = lift(slider_t.value) do n
    t[n]
end

for ax in axes[1:3]
    vlines!(ax, slice_pos)
end

slice1 = lift(slider_t.value) do n
    @view x1[n, :]
end
hist!(axes[4], slice1, normalization=:pdf, color=ColorSchemes.tab10[1])

slice2 = lift(slider_t.value) do n
    @view x2[n, :]
end
hist!(axes[5], slice2, normalization=:pdf, color=ColorSchemes.tab10[1])

slice3 = lift(slider_t.value) do n
    @view x3[n, :]
end
hist!(axes[6], slice3, normalization=:pdf, color=ColorSchemes.tab10[1])

xrange = -10:0.01:10
ylims!.(axes[4:6], 0, nothing)
xlims!.(axes[4:6], xrange[1], xrange[end])

pdf1 = lift(slider_t.value) do n
    normpdf(xrange, init_x, √(2 * D[1] * t[n]))
end
lines!(axes[4], xrange, pdf1, color=ColorSchemes.tab10[4])

pdf2 = lift(slider_t.value) do n
    normpdf(xrange, init_x, √(2 * D[2] * t[n]))
end
lines!(axes[5], xrange, pdf2, color=ColorSchemes.tab10[4])

pdf3 = lift(slider_t.value) do n
    normpdf(xrange, init_x, √(2 * D[3] * t[n]))
end
lines!(axes[6], xrange, pdf3, color=ColorSchemes.tab10[4])