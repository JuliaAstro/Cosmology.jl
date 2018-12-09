__precompile__()

module Cosmology

using QuadGK, Unitful
import Unitful: km, s
using UnitfulAstro: Mpc, Gpc, Gyr

export cosmology,
       age,
       angular_diameter_dist,
       comoving_radial_dist,
       comoving_transverse_dist,
       comoving_volume,
       comoving_volume_element,
       distmod,
       H,
       hubble_dist,
       hubble_time,
       luminosity_dist,
       lookback_time,
       scale_factor

abstract type AbstractCosmology end
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end

struct FlatLCDM{T<:Real} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_r))...)


a2E(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m*a + c.Ω_Λ*a^4)

struct ClosedLCDM{T<:Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    ClosedLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                       float(Ω_r))...)


struct OpenLCDM{T<:Real} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    OpenLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                     float(Ω_r))...)


function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a*a
    sqrt(c.Ω_r + c.Ω_m*a + (c.Ω_k + c.Ω_Λ*a2)*a2)
end

for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        struct $(name){T<:Real} <: $(Symbol("Abstract$(c)Cosmology"))
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

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    ade = exp((1 - 3*(c.w0 + c.wa))*log(a) + 3*c.wa*(a - 1))
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

    OmegaL = 1 - OmegaK - OmegaM - OmegaR

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
H(c::AbstractCosmology, z) = 100 * c.h * E(c, z) * km / s / Mpc

hubble_dist0(c::AbstractCosmology) = 2997.92458/c.h * Mpc
hubble_dist(c::AbstractCosmology, z) = hubble_dist0(c)/E(c,z)

hubble_time0(c::AbstractCosmology) = 9.77814/c.h * Gyr
hubble_time(c::AbstractCosmology, z) = hubble_time0(c)/E(c,z)

# distances

Z(c::AbstractCosmology, z::Real; kws...) =
    QuadGK.quadgk(a->1/a2E(c,a), scale_factor(z), 1; kws...)[1]

comoving_radial_dist(c::AbstractCosmology, z; kws...) = hubble_dist0(c)*Z(c, z; kws...)

comoving_transverse_dist(c::AbstractFlatCosmology, z; kws...) =
    comoving_radial_dist(c, z; kws...)
function comoving_transverse_dist(c::AbstractOpenCosmology, z; kws...)
    sqrtok = sqrt(c.Ω_k)
    hubble_dist0(c)*sinh(sqrtok*Z(c,z; kws...))/sqrtok
end
function comoving_transverse_dist(c::AbstractClosedCosmology, z; kws...)
    sqrtok = sqrt(abs(c.Ω_k))
    hubble_dist0(c)*sin(sqrtok*Z(c,z; kws...))/sqrtok
end

angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...)/(1 + z)

luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...)*(1 + z)

distmod(c::AbstractCosmology, z; kws...) =
    5 * log10(luminosity_dist(c, z; kws...) / Mpc) + 25

# volumes

comoving_volume(c::AbstractFlatCosmology, z; kws...) =
    (4pi/3)*(comoving_radial_dist(Gpc, c, z; kws...))^3
function comoving_volume(c::AbstractOpenCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...)/DH
    sqrtok = sqrt(c.Ω_k)
    2pi*(DH)^3*(x*sqrt(1 + c.Ω_k*x^2) - asinh(sqrtok*x)/sqrtok)/c.Ω_k
end
function comoving_volume(c::AbstractClosedCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c,z; kws...)/DH
    sqrtok = sqrt(abs(c.Ω_k))
    2pi*(DH)^3*(x*sqrt(1 + c.Ω_k*x^2) - asin(sqrtok*x)/sqrtok)/c.Ω_k
end

comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(Gpc, c)*angular_diameter_dist(Gpc, c, z; kws...)^2/a2E(c, scale_factor(z))

# times

T(c::AbstractCosmology, a0, a1; kws...) = QuadGK.quadgk(x->x/a2E(c,x), a0, a1; kws...)[1]
age(c::AbstractCosmology, z; kws...) = hubble_time0(c)*T(c, 0, scale_factor(z); kws...)
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c)*T(c, scale_factor(z), 1; kws...)

# Easily select a different unit
for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time)
    @eval $f(u::Unitful.Unitlike, args...; kws...) = uconvert(u, $f(args...; kws...))
end

###############
# Deprecations
#

@deprecate hubble_dist_mpc0(c::AbstractCosmology) ustrip(hubble_dist0(c::AbstractCosmology))
@deprecate hubble_dist_mpc(c::AbstractCosmology, z) ustrip(hubble_dist(c::AbstractCosmology, z))
@deprecate hubble_time_gyr0(c::AbstractCosmology) ustrip(hubble_time0(c::AbstractCosmology))
@deprecate hubble_time_gyr(c::AbstractCosmology, z) ustrip(hubble_time(c::AbstractCosmology, z))
@deprecate comoving_radial_dist_mpc(c::AbstractCosmology, z; kws...) ustrip(comoving_radial_dist(c::AbstractCosmology, z; kws...))
@deprecate comoving_transverse_dist_mpc(c::AbstractCosmology, z; kws...) ustrip(comoving_transverse_dist(c::AbstractCosmology, z; kws...))
@deprecate angular_diameter_dist_mpc(c::AbstractCosmology, z; kws...) ustrip(angular_diameter_dist(c::AbstractCosmology, z; kws...))
@deprecate luminosity_dist_mpc(c::AbstractCosmology, z; kws...) ustrip(luminosity_dist(c::AbstractCosmology, z; kws...))
@deprecate comoving_volume_gpc3(c::AbstractCosmology, z; kws...) ustrip(comoving_volume(c::AbstractCosmology, z; kws...))
@deprecate comoving_volume_element_gpc3(c::AbstractCosmology, z; kws...) ustrip(comoving_volume_element(c::AbstractCosmology, z; kws...))
@deprecate age_gyr(c::AbstractCosmology, z; kws...) ustrip(age(c::AbstractCosmology, z; kws...))
@deprecate lookback_time_gyr(c::AbstractCosmology, z; kws...) ustrip(lookback_time(c::AbstractCosmology, z; kws...))

end # module
