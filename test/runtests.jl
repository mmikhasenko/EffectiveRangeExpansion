using EffectiveRangeExpansion
using Parameters
using Test

import EffectiveRangeExpansion: fm_times_mev
#
const mπ⁰ = 0.1349768
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mD⁺ = 1.86965
# 
const mDˣ⁺ = 1.86483+145.4258e-3 # m(D) + Δm(D*,D) from PDG
const mDˣ⁰ = 2.00685
# 
const ΓDˣ⁺ = 83.4e-6
const ΓDˣ⁰ = 55.2e-6

e2m(e) = (mD⁰+mDˣ⁺)+e*1e-3
m2e(m) = (m-mD⁰-mDˣ⁺)*1e3

function k3b(e)
    m = e2m(e)
    M = sqrt(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺) # taken at the Dˣ⁺ pole
    p = cis(π/4)*sqrt((m-(M+mD⁰))*cis(-π/2))*
        sqrt(m+(M+mD⁰))*sqrt(m-(M-mD⁰))*sqrt(m+(M-mD⁰))/(2*m)  # branch cut down
    return p
end


@testset "Branch point in the complex plane" begin
    a₀_fm = -7.0 # fm
    r₀_fm =  3.0 # fm
    # 
    a⁻¹ = 1 / (a₀_fm / (1e-3*fm_times_mev))
    r = r₀_fm / (1e-3*fm_times_mev)
    #
    k = k3b
    D(e) = a⁻¹ + r/2*k(e)^2-1im*k(e)
    #
    # expansion is at 
    Eᵦ = m2e(sqrt(mDˣ⁺^2 - 1im * mDˣ⁺ * ΓDˣ⁺) + mD⁰)
    # 
    effrangepars = 
        effectiverangeexpansion(
            Δe->D(Δe+Eᵦ),
            Δe->k(Eᵦ+Δe),
            abs(imag(Eᵦ))/20)
    # 
    @unpack a_fm, r_fm = tophysicsunits(effrangepars)
    # 
    @test a_fm ≈ a₀_fm
    @test real(r_fm) ≈ r₀_fm
end