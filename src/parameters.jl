@with_kw struct ERP{T<:Number} # Effective Range Parameters
    N::T
    a⁻¹::T
    r::T
end
(f::ERP)(k) = f.N * (f.a⁻¹ + f.r * k^2 / 2 - 1im * k)
# 
function Base.show(io::IO, ere::ERP)
    println(io, typeof(ere), "(")
    println(io, "   a⁻¹= $(round(ere.a⁻¹, digits=3)),")
    println(io, "   r  = $(round(ere.r, digits=3)),")
    println(io, "   N  = $(round(ere.N, digits=3)))")
end
