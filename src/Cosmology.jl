module Cosmology

using QuadGK: quadgk
using Unitful
import Unitful: km, s, Gyr
using UnitfulAstro: Mpc, Gpc
using DocStringExtensions

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

"""
$(TYPEDEF)

Abstract supertype for all cosmological models.
"""
abstract type AbstractCosmology end
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end

"""
$(TYPEDEF)

ΛCDM model of the universe with ``Ω_k = 0``.
"""
struct FlatLCDM{T <: Real} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_r))...)

"""
$(TYPEDEF)

ΛCDM model of the universe with ``Ω_k < 0``.
"""
struct ClosedLCDM{T <: Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    ClosedLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                       float(Ω_r))...)

"""
$(TYPEDEF)

ΛCDM model of the universe with ``Ω_k > 0``.
"""
struct OpenLCDM{T <: Real} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    OpenLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                     float(Ω_r))...)


# define WCDM models, which includes a cosmological equation of state parameter w.
for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        struct $(name){T <: Real} <: $(Symbol("Abstract$(c)Cosmology"))
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

a2E(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m * a + c.Ω_Λ * a^4)
function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end
@doc raw"""
    a2E(c::Union{FlatLCDM,ClosedLCDM,OpenLCDM}, a)


Calculates the intermediate quantity ``a^2 E(a)``.
This is an internal function used to simplify computation.

Mathematical definition (for ΛCDM models):
```math
a^2 E(a) = \sqrt{Ω_r + Ω_m a + Ω_k a^2 + Ω_Λ a^4}
```
where ``Ω_k = 0`` for a flat cosmological model.
"""
a2E

@doc raw"""
    a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)

The implementation of ``a^2 E(a)`` for WCDM models.

Mathematical definition (for WCDM models):
```math
a^2 E(a) = \sqrt{Ω_r + Ω_m a + Ω_k a^2 + Ω_Λ a_{de}}
```
where ``a_{de} = \exp[(1 - 3 w_0 - 3 w_a) \log(a) + 3 w_a (a - 1)]``.
"""
function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    # dark energy scale factor
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end


"""
    cosmology(; h = 0.69,
                Neff = 3.04,
                OmegaK = 0,
                OmegaM = 0.29,
                OmegaR = nothing,
                Tcmb = 2.7255,
                w0 = -1,
                wa = 0)


# Parameters
* `h` - Dimensionless Hubble constant
* `Neff` - Effective number of massless neutrino species; used to compute Ω_ν
* `OmegaK` - Curvature density (Ω_k)
* `OmegaM` - Matter density (Ω_m)
* `OmegaR` - Radiation density (Ω_r)
* `Tcmb` - CMB temperature in Kelvin; used to compute Ω_γ
* `w0` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
* `wa` - CPL dark energy equation of state; `w = w0 + wa(1-a)`

# Examples
```jldoctest
julia> c = cosmology()
Cosmology.FlatLCDM{Float64}(0.69, 0.7099122024007928, 0.29, 8.77975992071536e-5)

julia> c = cosmology(OmegaK=0.1)
Cosmology.OpenLCDM{Float64}(0.69, 0.1, 0.6099122024007929, 0.29, 8.77975992071536e-5)

julia> c = cosmology(w0=-0.9, OmegaK=-0.1)
Cosmology.ClosedWCDM{Float64}(0.69, -0.1, 0.8099122024007929, 0.29, 8.77975992071536e-5, -0.9, 0.0)
```
"""
function cosmology(;h = 0.69,
                   Neff = 3.04,
                   OmegaK = 0,
                   OmegaM = 0.29,
                   OmegaR = nothing,
                   Tcmb = 2.7255,
                   w0 = -1,
                   wa = 0)

    if OmegaR === nothing
        OmegaG = 4.48131e-7 * Tcmb^4 / h^2
        OmegaN = Neff * OmegaG * (7 / 8) * (4 / 11)^(4 / 3)
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

"""
    scale_factor(z)

Return the scale factor ``a(t)`` for a given redshift ``z(t)``. According to the
[Friedmann–Lemaître–Robertson–Walker metric](https://en.wikipedia.org/wiki/Friedmann–Lemaître–Robertson–Walker_metric)
it's given as ``a = 1/(1 + z)`` ([Schneider 2015, p. 186](@cite schneider2015)).

A scale factor of 1, i.e., a redshift of 0, refers to the present epoch.
"""
scale_factor(z) = 1 / (1 + z)

@doc raw"""
    E(c::AbstractCosmology, z)

Dimensionless Hubble function ``E(z)`` at redshift `z`. It's defined as
```math
E(z) ≡ \frac{H(z)}{H_0} = \frac{H(z)}{(100\mathrm{km/s/Mpc}) h}
```
where ``H_0 = H(z=0)`` is the Hubble parameter at the present epoch
([Schneider 2015, p. 183](@cite schneider2015)).
"""
E(c::AbstractCosmology, z) = (a = scale_factor(z); a2E(c, a) / a^2)

"""
    H(c::AbstractCosmology, z)

Hubble parameter at redshift `z`.
"""
H(c::AbstractCosmology, z) = 100 * c.h * E(c, z) * km / s / Mpc

"""
    hubble_dist0(c::AbstractCosmology)

Hubble distance at redshift 0.

### See also
[`hubble_dist`](@ref)
"""
hubble_dist0(c::AbstractCosmology) = 2997.92458 / c.h * Mpc
"""
    hubble_dist(c::AbstractCosmology, z)

Hubble distance ``D_H``, defined as the product of the speed of light and
the Hubble time. That is, ``D_H(z) = c / H(z)``.

### See also
[`hubble_time`](@ref)
"""
hubble_dist(c::AbstractCosmology, z) = hubble_dist0(c) / E(c, z)

"""
    hubble_time0(c::AbstractCosmology)

Hubble time at redshift 0.

### See also
[`hubble_time`](@ref)
"""
hubble_time0(c::AbstractCosmology) = 9.777922216807891 / c.h * Gyr
"""
    hubble_time(c::AbstractCosmology, z)

Hubble time, defined as the inverse of the Hubble parameter. That is,
``t_H(z) = 1/H(z)``.

### See also
[`hubble_dist`](@ref)
"""
hubble_time(c::AbstractCosmology, z) = hubble_time0(c) / E(c, z)

# distances

Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) =
    quadgk(a->1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) =
    quadgk(a->1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]
@doc raw"""
    Z(c::AbstractCosmology, z, nothing; kws...)
    Z(c::AbstractCosmology, z₁, z₂; kws...)

The line-of-sight comoving distance contributions for comoving radial distance.

It performs the integral
```math
Z = \int_{z_1}^{z_2} \frac{dz}{E(z)} = \int_{a_2}^{a_1} \frac{da}{a^2 E(a)}
```
where we can perform a change of variables with ``a = 1/(1+z)``,
and ``dz = -da/a^2``.

If `nothing` is used for the second bound of integration, it defaults
to `z₁ = 0` (i.e., `a₁ = 1`).

### See also
[`comoving_radial_dist`](@ref)
"""
Z

comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

@doc raw"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Comoving radial distance (``D_C``) in Mpc at redshift `z₂` as seen by an observer at `z₁`.
Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if
provided.

It's calculated as ``D_C = D_{H0} Z``, where ``D_{H0}`` is the Hubble distance at
the present epoch and, ``Z = \int_{z_1}^{z_2} \frac{dz}{E(z)}``.
"""
comoving_radial_dist


comoving_transverse_dist(c::AbstractFlatCosmology, z₁, z₂ = nothing; kws...) =
    comoving_radial_dist(c, z₁, z₂; kws...)
function comoving_transverse_dist(c::AbstractOpenCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(c.Ω_k)
    hubble_dist0(c) * sinh(sqrtok * Z(c, z₁, z₂; kws...)) / sqrtok
end
function comoving_transverse_dist(c::AbstractClosedCosmology, z₁, z₂ = nothing; kws...)
    sqrtok = sqrt(abs(c.Ω_k))
    hubble_dist0(c) * sin(sqrtok * Z(c, z₁, z₂; kws...)) / sqrtok
end

angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) / (1 + z)
angular_diameter_dist(c::AbstractCosmology, z₁, z₂; kws...) =
    comoving_transverse_dist(c, z₁, z₂; kws...) / (1 + z₂)

"""
    angular_diameter_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Ratio of the proper transverse size in Mpc of an object at redshift `z₂` to its
angular size in radians, as seen by an observer at `z₁`.  Redshift `z₁` defaults
to 0 if omitted.  Will convert to compatible unit `u` if provided.
"""
angular_diameter_dist

luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to
compatible unit `u` if provided.
"""
luminosity_dist

"""
    distmod(c::AbstractCosmology, z)

Distance modulus in magnitudes at redshift `z`.
"""
distmod(c::AbstractCosmology, z; kws...) =
    5 * log10(luminosity_dist(c, z; kws...) / Mpc) + 25

# volumes

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z)

Comoving volume in cubic Gpc out to redshift `z`. Will convert to compatible unit `u` if provided.
"""
comoving_volume(c::AbstractFlatCosmology, z; kws...) =
    (4pi / 3) * (comoving_radial_dist(Gpc, c, z; kws...))^3
function comoving_volume(c::AbstractOpenCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...) / DH
    sqrtok = sqrt(c.Ω_k)
    2pi * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asinh(sqrtok * x) / sqrtok) / c.Ω_k
end
function comoving_volume(c::AbstractClosedCosmology, z; kws...)
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...) / DH
    sqrtok = sqrt(abs(c.Ω_k))
    2pi * (DH)^3 * (x * sqrt(1 + c.Ω_k * x^2) - asin(sqrtok * x) / sqrtok) / c.Ω_k
end

"""
    comoving_volume_element([u::Unitlike,] c::AbstractCosmology, z)

Comoving volume element in Gpc out to redshift `z`. Will convert to compatible unit `u` if provided.
"""
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(Gpc, c) * angular_diameter_dist(Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

# times

T(c::AbstractCosmology, a0, a1; kws...) = quadgk(x->x / a2E(c, x), a0, a1; kws...)[1]

"""
    age([u::Unitlike,] c::AbstractCosmology, z)

Age of the universe in Gyr at redshift `z`. Will convert to compatible unit `u` if provided.
"""
age(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, 0, scale_factor(z); kws...)

"""
    lookback_time([u::Unitlike,] c::AbstractCosmology, z)

Difference between age at redshift 0 and age at redshift `z` in Gyr.
Will convert to compatible unit `u` if provided.
"""
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, scale_factor(z), 1; kws...)

# Easily select a different unit
for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time,
          :comoving_radial_dist, :comoving_transverse_dist,
          :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element,
          :age, :lookback_time)
    @eval $f(u::Unitful.Unitlike, args...; kws...) = uconvert(u, $f(args...; kws...))
end

end # module
