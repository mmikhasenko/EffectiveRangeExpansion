# Effective Range Expansion

The package provides tools to numerically estimate the scattering length and effective range for given amplitude.

Details of the calculation are described in [this paper](https://inspirehep.net/literature/2048996).

## Definition
The expansion series reads:

$$f(s) = N\bigg( a^{-1} + r \frac{k^2}{2} - i k(s) + o(k^4)\bigg)$$




where `f` is the target function, `N` is a numerical constant, `k` is the break-up momentum in the system that vanishes at the expansion point (the branch point).
The `a` and `r` are the scattering length and the effective range.

## Computation
The computation is done by computing cricular (cauchy) integrals about the expansion point and matching their result for the
left and the right side of the Eq. (1).

Use the method,
```julia
effectiverangeexpansion(f, k, method::EffectiveRangeExpansionMethod)
```
to obtained a named tuple `(; a⁻¹, r, N)`.

The method parameter adjusts the algorithm. The `ComplexBranchPointExpansion(CircularIntegral(r))` is the default option.
A discrete version of the integral is available to speed up the calculation.
```julia
struct CircularIntegral <: CircularIntegralMethod
    r::Real
end
#
struct CircularSum <: CircularIntegralMethod
    r::Real
    N::Int
end
```

