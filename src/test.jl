using Cosmology
using Base.Test

function test_approx_eq_rtol(va, vb, rtol, astr, bstr)
    diff = max(abs(va - vb))
    tol = rtol*max(max(abs(va)), max(abs(vb)))
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
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6607.6579 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.4694 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.7527 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.4694-5.7527 age_rtol

c = cosmology(h=0.7, OmegaK=0.1, OmegaM=0.3, OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1619.9588 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3209.784 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6479.8352 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.064 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.5466 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.064-5.5466 age_rtol

c = cosmology(h=0.7, OmegaK=-0.1, OmegaM=0.3, OmegaR=0)
@test_approx_eq_rtol angular_diameter_dist_mpc(c,1) 1686.5272 dist_rtol
@test_approx_eq_rtol comoving_radial_dist_mpc(c,1) 3408.937 dist_rtol
@test_approx_eq_rtol luminosity_dist_mpc(c,1) 6746.1088 dist_rtol
@test_approx_eq_rtol age_gyr(c,0) 13.925 age_rtol
@test_approx_eq_rtol age_gyr(c,1) 5.9868 age_rtol
@test_approx_eq_rtol lookback_time_gyr(c,1) 13.925-5.9868 age_rtol
