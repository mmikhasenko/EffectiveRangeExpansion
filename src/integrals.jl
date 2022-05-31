circleintegral(f, x₀::Real, r::Real) =
    quadgk(ϕ -> f(x₀ + r * cis(ϕ)) * r * cis(ϕ), -7π / 8, 9π / 8)[1] / (2π) # cancel 1im
# 
circleintegral(f, x₀::Real, r::Real, N::Int) =
    sum(f(x₀ + r * cis(ϕ)) * r * cis(ϕ) for ϕ in range(-7π / 8, 9π / 8, length=N + 1)[2:end])[1] / N # 2π i cancels
#
cauchyintegral(F, x₀::Real, r) = circleintegral(x′ -> F(x′ + x₀) / x′, 0.0, r)
cauchyintegral′(F, x₀::Real, r) = circleintegral(x′ -> F(x′ + x₀) / x′^2, 0.0, r)
cauchyintegral′′(F, x₀::Real, r) = circleintegral(x′ -> 2F(x′ + x₀) / x′^3, 0.0, r)
# 


abstract type CircularIntegralMethod end
struct CircularIntegral <: CircularIntegralMethod
    r::Real
end
struct CircularSum <: CircularIntegralMethod
    r::Real
    N::Int
end

circleintegral(F, x₀::Real, method::CircularIntegral) = circleintegral(F, x₀, method.r)
circleintegral(F, x₀::Real, method::CircularSum) = circleintegral(F, x₀, method.r, method.N)

