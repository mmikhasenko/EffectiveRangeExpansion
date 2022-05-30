
const fm_times_mev = 197.3269804

function tophysicsunits(p::NamedTuple)
    @unpack a⁻¹, r = p
    a_fm = 1e-3 * fm_times_mev / real(p.a⁻¹)
    r_fm = 1e-3 * fm_times_mev * p.r
    (; a_fm, r_fm)
end