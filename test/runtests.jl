using Cosmology
using Test, Unitful, UnitfulAstro, QuadGK

# values from http://icosmos.co.uk/

dist_rtol = 1e-6
age_rtol = 2e-4
# Integrating a unitful function would require UnitfulIntegration.jl.  Without using it, we
# strip the units away from the integrand function
integrand(c, z) = 4pi*ustrip(comoving_volume_element(c, z))

@testset "FlatLCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9145u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3303.829u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 151.0571u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6607.6579u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1002 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.4694u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.7527u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.4694-5.7527)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "OpenLCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1619.9588u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3209.784u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.0856u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6479.8352u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0578 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.064u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.5466u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.064-5.5466)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "ClosedLCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1686.5272u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3408.937u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 163.8479u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6746.1088u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1453 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.925u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.9868u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.925-5.9868)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "FlatWCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.0585u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3224.1169u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.3851u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6448.2338u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0472 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.1915u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.6464u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.1915-5.6464)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "OpenWCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1588.0181u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,rtol=dist_rtol,1) ≈ 3147.6227u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 132.0466u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6352.0723u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0146 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 12.8488u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.4659u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (12.8488-5.4659)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "ClosedWCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1637.5993u"Mpc" rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3307.9932u"Mpc" rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 149.8301u"Gpc^3" rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ ustrip(comoving_volume(c, 2.5))
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6550.3973u"Mpc" rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0813 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.5702u"Gyr" rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.8482u"Gyr" rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.5702-5.8482)u"Gyr" rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "Non-Float64" begin
    # Test that FlatLCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=0.7, OmegaM=big(0.3), OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9145u"Mpc" rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.4030879e10u"Mpc^3" rtol = dist_rtol
    # Test that FlatWCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=big(0.7), OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.0585u"Mpc" rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.1378625e10u"Mpc^3" rtol = dist_rtol
end

@testset "Unit conversion" begin
    c = cosmology(h=0.9, OmegaM=0.5, OmegaR=0)
    for u in (u"m", u"pc", u"ly")
        @test unit(luminosity_dist(u, c, 1)) == u
        @test unit(angular_diameter_dist(u, c, 2)) == u
    end
    for u in (u"s", u"yr")
        @test unit(age(u, c, 3)) == u
        @test unit(lookback_time(u, c, 4)) == u
    end
end

@testset "Utilities" begin
    c = cosmology(h = 0.7)
    @test hubble_time(c, 0) ≈ Cosmology.hubble_time0(c)
    @test hubble_dist(c, 0) ≈ Cosmology.hubble_dist0(c)
    @test H(c, 0) ≈ 70u"km/s/Mpc"
end
