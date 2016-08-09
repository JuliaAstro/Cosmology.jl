module Cosmology

using Compat

export cosmology,
       age_gyr,
       angular_diameter_dist_mpc,
       comoving_radial_dist_mpc,
       comoving_transverse_dist_mpc,
       comoving_volume_gpc3,
       distmod,
       H,
       hubble_dist_mpc,
       hubble_time_gyr,
       luminosity_dist_mpc,
       lookback_time_gyr,
       scale_factor

abstract AbstractCosmology
abstract AbstractClosedCosmology <: AbstractCosmology
abstract AbstractFlatCosmology <: AbstractCosmology
abstract AbstractOpenCosmology <: AbstractCosmology

immutable FlatLCDM{T<:Real} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_r))...)


a2E(c::FlatLCDM, a::Float64) = sqrt(c.Ω_r + c.Ω_m*a + c.Ω_Λ*a^4)

immutable ClosedLCDM{T<:Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    ClosedLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                       float(Ω_r))...)


immutable OpenLCDM{T<:Real} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    OpenLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                     float(Ω_r))...)


function a2E(c::Union{ClosedLCDM,OpenLCDM}, a::Float64)
    a2 = a*a
    sqrt(c.Ω_r + c.Ω_m*a + (c.Ω_k + c.Ω_Λ*a2)*a2)
end

for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        immutable $(name){T<:Real} <: $(Symbol("Abstract$(c)Cosmology"))
            h::T
            Ω_k::T
            Ω_Λ::T
            Ω_m::T
            Ω_r::T
            w0::T
            wa::T
        end
        function $(name)(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real,
                         w0::Real, wa::Real)
            $(name)(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                            float(Ω_r), float(w0), float(wa))...)
        end
    end
end

function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real, w0::Real, wa::Real)
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    end
end

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a::Float64)
    ade = exp((1.0 - 3.0*(c.w0 + c.wa))*log(a) + 3.0*c.wa*(a - 1.0))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k*a)*a + c.Ω_Λ*ade)
end

function cosmology(;h=0.69,
                   Neff=3.04,
                   OmegaK=0,
                   OmegaM=0.29,
                   OmegaR=nothing,
                   Tcmb=2.7255,
                   w0=-1,
                   wa=0)

    if OmegaR === nothing
        OmegaG = 4.48131e-7*Tcmb^4/h^2
        OmegaN = Neff*OmegaG*(7/8)*(4/11)^(4/3)
        OmegaR = OmegaG + OmegaN
    end

    OmegaL = 1. - OmegaK - OmegaM - OmegaR

    if !(w0 == -1 && wa == 0)
        return WCDM(h, OmegaK, OmegaL, OmegaM, OmegaR, w0, wa)
    end

    if OmegaK < 0
        return ClosedLCDM(h, OmegaK, OmegaL, OmegaM, OmegaR)
    elseif OmegaK > 0
        return OpenLCDM(h, OmegaK, OmegaL, OmegaM, OmegaR)
    else
        return FlatLCDM(h, OmegaL, OmegaM, OmegaR)
    end
end

# hubble rate

scale_factor(z) = 1/(1 + z)
E(c::AbstractCosmology, z) = (a = scale_factor(z); a2E(c,a)/a^2)
H(c::AbstractCosmology, z) = 100. * c.h * E(c, z)

hubble_dist_mpc0(c::AbstractCosmology) = 2997.92458/c.h
hubble_dist_mpc(c::AbstractCosmology, z) = hubble_dist_mpc0(c)/E(c,z)

hubble_time_gyr0(c::AbstractCosmology) = 9.77814/c.h
hubble_time_gyr(c::AbstractCosmology, z) = hubble_time_gyr0(c)/E(c,z)

# distances

Z(c::AbstractCosmology, z::Real) = ((q,_) = quadgk(a::Float64->1.0/a2E(c,a), scale_factor(z), 1); q)

comoving_radial_dist_mpc(c::AbstractCosmology, z) = hubble_dist_mpc0(c)*Z(c, z)

comoving_transverse_dist_mpc(c::AbstractFlatCosmology, z) =
    comoving_radial_dist_mpc(c, z)
function comoving_transverse_dist_mpc(c::AbstractOpenCosmology, z)
    sqrtok = sqrt(c.Ω_k)
    hubble_dist_mpc0(c)*sinh(sqrtok*Z(c,z))/sqrtok
end
function comoving_transverse_dist_mpc(c::AbstractClosedCosmology, z)
    sqrtok = sqrt(abs(c.Ω_k))
    hubble_dist_mpc0(c)*sin(sqrtok*Z(c,z))/sqrtok
end

angular_diameter_dist_mpc(c::AbstractCosmology, z) =
    comoving_transverse_dist_mpc(c, z)/(1 + z)

luminosity_dist_mpc(c::AbstractCosmology, z) =
    comoving_transverse_dist_mpc(c, z)*(1 + z)

distmod(c::AbstractCosmology, z) =
    5.0 * log10(luminosity_dist_mpc(c, z)) + 25.0

# volumes

comoving_volume_gpc3(c::AbstractFlatCosmology, z) =
    (4pi/3)*(comoving_radial_dist_mpc(c,z)*1e-3)^3
function comoving_volume_gpc3(c::AbstractOpenCosmology, z)
    DH = hubble_dist_mpc0(c)
    x = comoving_transverse_dist_mpc(c,z)/DH
    sqrtok = sqrt(c.Ω_k)
    2pi*(DH*1e-3)^3*(x*sqrt(1. + c.Ω_k*x^2) - asinh(sqrtok*x)/sqrtok)/c.Ω_k
end
function comoving_volume_gpc3(c::AbstractClosedCosmology, z)
    DH = hubble_dist_mpc0(c)
    x = comoving_transverse_dist_mpc(c,z)/DH
    sqrtok = sqrt(abs(c.Ω_k))
    2pi*(DH*1e-3)^3*(x*sqrt(1. + c.Ω_k*x^2) - asin(sqrtok*x)/sqrtok)/c.Ω_k
end

comoving_volume_element_gpc3(c::AbstractCosmology, z) =
    1e-9*hubble_dist_mpc0(c,z)*angular_diameter_dist_mpc(c,z)^2/a2E(c,scale_factor(z))

# times

T(c::AbstractCosmology, a0::Float64, a1::Float64) = ((q,_) = quadgk(x::Float64->x/a2E(c,x), a0, a1); q)
age_gyr(c::AbstractCosmology, z) = hubble_time_gyr0(c)*T(c, 0., scale_factor(z))
lookback_time_gyr(c::AbstractCosmology, z) = hubble_time_gyr0(c)*T(c, scale_factor(z), 1.)

end # module
