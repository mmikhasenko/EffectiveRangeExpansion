
abstract type EffectiveRangeExpansionMethod end
struct ComplexBranchPointExpansion <: EffectiveRangeExpansionMethod
    cim::CircularIntegralMethod
end

effectiverangeexpansion(f, k, r::Float64) =
    effectiverangeexpansion(f, k, ComplexBranchPointExpansion(CircularIntegral(r)))

"""
f(s) / N = a⁻¹ + r k(s)^2 / 2 - i k(s)

Is calculated using cauchy integral theorem.
"""
function effectiverangeexpansion(f, k, method::EffectiveRangeExpansionMethod)
    N = circleintegral(f, 0.0, method.cim) /
        circleintegral(k, 0.0, method.cim) / (-1im)
    # 
    f̂₀ = cauchyintegral(x -> f(x) / N + 1im * k(x), 0.0, method.cim)
    f̂′ = cauchyintegral′(x -> f(x) / N + 1im * k(x), 0.0, method.cim)
    # 
    # k²₀  = cauchyintegral(  x->k(x)^2, 0.0, r)
    k²′ = cauchyintegral′(x -> k(x)^2, 0.0, method.cim)
    # 
    a⁻¹ = f̂₀
    r = 2 * f̂′ / k²′
    #
    return (; a⁻¹, r, N)
end