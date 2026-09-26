using Test
using Cosmology
using QuadGK
using DynamicQuantities: @u_str, @us_str, dimension, ustrip
using DynamicQuantities.Constants: Mpc, Gpc
using DynamicQuantities.Units: Gyr

# values from http://icosmos.co.uk/

dist_rtol = 1e-6
age_rtol = 2e-4
integrand(c, z) = 4π * comoving_volume_element(c, z)

@testset "FlatLCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9145 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 625.3444 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3303.829 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 151.0571 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6607.6579 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1002 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.4694 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.7527 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.4694-5.7527) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "OpenLCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1619.9588 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 598.9118 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3209.784 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.0856 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6479.8352 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0578 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.064 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.5466 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.064-5.5466) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "ClosedLCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1686.5272 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 655.6019 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3408.937 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 163.8479 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6746.1088 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1453 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.925 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.9868 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.925-5.9868) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "FlatWCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.0585 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 607.6802 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3224.1169 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 140.3851 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6448.2338 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0472 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.1915 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.6464 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.1915-5.6464) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "OpenWCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1588.0181 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 585.4929 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,rtol=dist_rtol,1) ≈ 3147.6227 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 132.0466 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6352.0723 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0146 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 12.8488 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.4659 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (12.8488-5.4659) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "ClosedWCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1637.5993 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 632.5829 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_radial_dist(c,1,rtol=dist_rtol) ≈ 3307.9932 * Mpc rtol = dist_rtol
    @test comoving_volume(c,1,rtol=dist_rtol) ≈ 149.8301 * Gpc^3 rtol = dist_rtol
    @test quadgk(z -> integrand(c, z), 0, 2.5)[1] ≈ comoving_volume(c, 2.5)
    @test luminosity_dist(c,1,rtol=dist_rtol) ≈ 6550.3973 * Mpc rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0813 rtol = dist_rtol
    @test age(c,0,rtol=age_rtol) ≈ 13.5702 * Gyr rtol = age_rtol
    @test age(c,1,rtol=age_rtol) ≈ 5.8482 * Gyr rtol = age_rtol
    @test lookback_time(c,1,rtol=age_rtol) ≈ (13.5702-5.8482) * Gyr rtol = age_rtol
    @test age(c, 1) + lookback_time(c, 1) ≈ age(c, 0)
end

@testset "Non-Float64" begin
    # Test that FlatLCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=0.7, OmegaM=big(0.3), OmegaR=0)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1651.9145 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 625.3444 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.4030879e10 * Mpc^3 rtol = dist_rtol
    # Test that FlatWCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=big(0.7), OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist(c,1,rtol=dist_rtol) ≈ 1612.0585 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,1,2,rtol=dist_rtol) ≈ 607.6802 * Mpc rtol = dist_rtol
    @test angular_diameter_dist(c,pi,rtol=dist_rtol) ≈ angular_diameter_dist(c,0,pi,rtol=dist_rtol) rtol = dist_rtol
    @test comoving_volume_element(c, big(1.41)) ≈ 3.1378625e10 * Mpc^3 rtol = dist_rtol
end

@testset "Unit conversion" begin
    c = cosmology(h=0.9, OmegaM=0.5, OmegaR=0)
    for u in (us"m", us"Constants.pc", us"Constants.ly")
        @test dimension(luminosity_dist(u, c, 1)) == dimension(u)
        @test dimension(angular_diameter_dist(u, c, 2)) == dimension(u)
    end
    for u in (us"s", us"yr")
        @test dimension(age(u, c, 3)) == dimension(u)
        @test dimension(lookback_time(u, c, 4)) == dimension(u)
    end
end
