# Using QuadGK, we can numerically integrate the Gaussian integrals in the assignment
#
using QuadGK
# Define the integrand
μ = 1.0
σ = 1.0
f(x) = exp(-(x-μ)^2/(2σ^2))

println(quadgk(f, -Inf, Inf, rtol=1e-9))
println(quadgk(f, 0, Inf, rtol=1e-9))

