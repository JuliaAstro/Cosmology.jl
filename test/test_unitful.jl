using Cosmology
using Test, Unitful, UnitfulAstro
using DynamicQuantities: DynamicQuantities as DQ
using DynamicQuantities.Constants: Mpc, Gpc
using DynamicQuantities.Units: Gyr


@testset "interface"  begin
    cosmo = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)
    for (f, u, u_DQ) in (
        (:age, u"Gyr", Gyr),
        (:angular_diameter_dist, u"Mpc", Mpc),
        (:comoving_radial_dist, u"Mpc", Mpc),
        (:comoving_transverse_dist, u"Mpc", Mpc),
        (:comoving_volume, u"Gpc^3", Gpc^3),
        (:comoving_volume_element, u"Gpc^3", Gpc^3),
        (:hubble_dist, u"Mpc", Mpc),
        (:hubble_time, u"Gyr", Gyr),
        (:luminosity_dist, u"Mpc", Mpc),
        (:lookback_time, u"Gyr", Gyr),
        )
        @eval begin
            q = $f($cosmo, 1.0)
            q_unitful = $f($u, $cosmo, 1.0)
            @test eltype(q_unitful) <: Unitful.Quantity
            @test ustrip($u, q_unitful) ≈ DQ.ustrip($u_DQ, q)
        end
    end
end
