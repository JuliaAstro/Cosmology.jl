include("../src/Cosmology.jl")
using Cosmology
using Base.Test

function test_approx_eq_rtol(va, vb, rtol, astr, bstr)
    diff = maximum(abs(va - vb))
    tol = rtol*max(maximum(abs(va)), maximum(abs(vb)))
    if diff > tol
        sdiff = string("|", astr, " - ", bstr, "| <= ", tol)
        error("assertion failed: ", sdiff,
              "\n  ", astr, " = ", va,
              "\n  ", bstr, " = ", vb,
              "\n  difference = ", diff, " > ", tol)
    end
end

macro test_approx_eq_rtol(a, b, c)
    :(test_approx_eq_rtol($(esc(a)), $(esc(b)), $(esc(c)), $(string(a)), $(string(b))))
end

# values from http://icosmos.co.uk/

dist_rtol = 1e-6
age_rtol = 2e-4

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1651.9145 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3303.829 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 151.0571 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6607.6579 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.1002 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.4694 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.7527 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.4694-5.7527 age_rtol

c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1619.9588 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3209.784 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 140.0856 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6479.8352 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.0578 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.064 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.5466 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.064-5.5466 age_rtol

c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1686.5272 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3408.937 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 163.8479 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6746.1088 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.1453 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.925 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.9868 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.925-5.9868 age_rtol

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1612.0585 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3224.1169 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 140.3851 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6448.2338 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.0472 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.1915 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.6464 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.1915-5.6464 age_rtol

c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1588.0181 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3147.6227 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 132.0466 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6352.0723 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.0146 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 12.8488 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.4659 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 12.8488-5.4659 age_rtol

c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1637.5993 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3307.9932 dist_rtol
@test_approx_eq_rtol comoving_volume_gpc3(c,1) 149.8301 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6550.3973 dist_rtol
@test_approx_eq_rtol distmod(c,1) 44.0813 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.5702 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.8482 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.5702-5.8482 age_rtol

# Test that FlatLCDM works with non-Float64 (BigFloat in this example)
c = cosmology(h=0.7, OmegaM=big(0.3), OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1651.9145 dist_rtol

# Test that FlatWCDM works with non-Float64 (BigFloat in this example)
c = cosmology(h=big(0.7), OmegaM=0.3, OmegaR=0, w0=-0.9, wa=0.1)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1612.0585 dist_rtol
