using Cosmology
using Test

# values from http://icosmos.co.uk/

dist_rtol = 1e-6
age_rtol = 2e-4

@testset "FlatLCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1651.9145 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3303.829 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 151.0571 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6607.6579 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1002 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 13.4694 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.7527 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 13.4694-5.7527 rtol = age_rtol
end

@testset "OpenLCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1619.9588 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3209.784 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 140.0856 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6479.8352 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0578 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 13.064 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.5466 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 13.064-5.5466 rtol = age_rtol
end

@testset "ClosedLCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1686.5272 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3408.937 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 163.8479 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6746.1088 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.1453 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 13.925 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.9868 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 13.925-5.9868 rtol = age_rtol
end

@testset "FlatWCDM" begin
    c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1612.0585 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3224.1169 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 140.3851 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6448.2338 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0472 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 13.1915 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.6464 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 13.1915-5.6464 rtol = age_rtol
end

@testset "OpenWCDM" begin
    c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1588.0181 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3147.6227 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 132.0466 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6352.0723 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0146 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 12.8488 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.4659 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 12.8488-5.4659 rtol = age_rtol
end

@testset "ClosedWCDM" begin
    c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1637.5993 rtol = dist_rtol
    @test comoving_radial_dist_mpc(c,1,rtol=dist_rtol) ≈ 3307.9932 rtol = dist_rtol
    @test comoving_volume_gpc3(c,1,rtol=dist_rtol) ≈ 149.8301 rtol = dist_rtol
    @test luminosity_dist_mpc(c,1,rtol=dist_rtol) ≈ 6550.3973 rtol = dist_rtol
    @test distmod(c,1,rtol=dist_rtol) ≈ 44.0813 rtol = dist_rtol
    @test age_gyr(c,0,rtol=age_rtol) ≈ 13.5702 rtol = age_rtol
    @test age_gyr(c,1,rtol=age_rtol) ≈ 5.8482 rtol = age_rtol
    @test lookback_time_gyr(c,1,rtol=age_rtol) ≈ 13.5702-5.8482 rtol = age_rtol
end

@testset "Non-Float64" begin
    # Test that FlatLCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=0.7, OmegaM=big(0.3), OmegaR=0)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1651.9145 rtol = dist_rtol
    # Test that FlatWCDM works with non-Float64 (BigFloat in this example)
    c = cosmology(h=big(0.7), OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
    @test angular_diameter_dist_mpc(c,1,rtol=dist_rtol) ≈ 1612.0585 rtol = dist_rtol
end
