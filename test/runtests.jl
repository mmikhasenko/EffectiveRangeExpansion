using EffectiveRangeExpansion
using Parameters
using Test

import EffectiveRangeExpansion: fm_times_mev


@testset "Break-up Rotate Cut" begin
    a,b = breakuprotatecut.((3-0.1im) .+ 1e-3 .*[-1,1]; m1=1, m2=2, ϕ=-π / 2)
    @test isapprox(-a, b; rtol=0.02) # 2%
end



@testset "ERP implementes a⁻¹ + r/2*k^2-1im*k" begin
    a⁻¹, r = 0.1, 4.2
    erp0 = ERP(; a⁻¹, r, N=1.0) # a⁻¹ + r/2*k^2-1im*k
    @test erp0(0) == a⁻¹
    @test erp0(1)-erp0(-1) == -2im
end


@testset "Branch point in the complex plane" begin
    a₀_fm = -7.0 # fm
    r₀_fm =  3.0 # fm
    # 
    a⁻¹ = 1 / (a₀_fm / (1e-3*fm_times_mev))
    r = r₀_fm / (1e-3*fm_times_mev)
    #
    mD⁰ = 1.86483
    mDˣ⁺ = 1.86483+145.4258e-3 # m(D) + Δm(D*,D) from PDG
    ΓDˣ⁺ = 83.4e-6

    e2m(e) = (mD⁰+mDˣ⁺)+e*1e-3
    k(e) = breakuprotatecut(e2m(e);
        m1=sqrt(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺), m2=mD⁰)
    # 
    D(e) = a⁻¹ + r/2*k(e)^2-1im*k(e)
    #
    # expansion is at
    m2e(m) = (m - mD⁰-mDˣ⁺)*1e3
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