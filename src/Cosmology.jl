module Cosmology

using QuadGK: quadgk
using Unitful
import Unitful: km, s, Gyr
using UnitfulAstro: Mpc, Gpc
using OrdinaryDiffEq
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
       scale_factor,
       f_DE,
       growth_factor

"""
$(TYPEDEF)

Abstract supertype for all cosmological models.
"""
abstract type AbstractCosmology end

# define all the ΛCDM and wCDM models, the latter of which includes a
# cosmological equation of state parameter w.
for (model, prettyname) in (("LCDM", "ΛCDM"), ("WCDM", "wCDM"))
    @eval begin
        """
        $(TYPEDEF)

        $($prettyname) model of the universe.
        """
        Base.@kwdef struct $model{T <: Real} <: AbstractCosmology
            h::T = 0.67
            Ω_k::T = 0.0
            Ω_c::T = 0.3
            Ω_b::T = 0.0
            Neff::T = 3.04
            Tcmb::T = 2.755
            w0::T = -1.0
            wa::T = 0.0

            # Derived densities
            Ω_γ::T = 4.48131e-7 * Tcmb^4 / h^2
            Ω_ν::T = Neff * Ω_γ * (7 / 8) * (4 / 11)^(4 / 3)
            Ω_r::T = Ω_γ + Ω_ν
            Ω_m::T = Ω_c + Ω_b
            Ω_Λ::T = 1.0 - Ω_m - Ω_r - Ω_k
        end
        function $model(h, Ω_k, Ω_c, Ω_b, Neff, Tcmb, w0, wa, Ω_γ, Ω_ν, Ω_r, Ω_m, Ω_Λ)
            return $model(promote(h, Ω_k, Ω_c, Ω_b, Neff, Tcmb, w0, wa, Ω_γ, Ω_ν, Ω_r, Ω_m, Ω_Λ)...)
        end
        #$model(; kwargs...) = $model(; promote(map(float, kwargs)...)...)
    end
end
"""
$(TYPEDEF)

ΛCDM model of the universe.
"""
struct LCDM{T <: Real} <: AbstractCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
LCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    LCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m), float(Ω_r))...)

"""
$(TYPEDEF)

wCDM model of the universe, which includes a cosmological equation of state parameter w.
"""
struct WCDM{T <: Real} <: AbstractCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
    w0::T
    wa::T
end
function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real, w0::Real, wa::Real)
    return WCDM(
        promote(
            float(h), float(Ω_k), float(Ω_Λ), float(Ω_m), float(Ω_r),
            float(w0), float(wa),
        )...
    )
end

@doc raw"""
    a2E(c::AbstractCosmology, a)

Calculates the intermediate quantity ``a^2 E(a)``.
This is an internal function used to simplify computation.

Mathematical definition for ΛCDM models:
```math
a^2 E(a) = \sqrt{Ω_r + Ω_m a + Ω_k a^2 + Ω_Λ a^4}
```
and for wCDM models:
```math
a^2 E(a) = \sqrt{Ω_r + Ω_m a + Ω_k a^2 + Ω_Λ a_{de}}
```
where ``Ω_k = 0`` for a flat cosmological model,
and ``a_{de} = a^{1 - 3(w_0 + w_a)} \exp(3 w_a (a - 1))`` [Scherrer2015](@cite).
"""
function a2E end
function a2E(c::LCDM, a)
    a2 = a * a
    return sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end
a2E(c::WCDM, a) = sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade(c, a))

# dark energy scale factor
ade(c::WCDM, a) = a^(1 - 3 * (c.w0 + c.wa)) * exp(3 * c.wa * (a - 1))

f_DE(c::Union{FlatLCDM, ClosedLCDM, OpenLCDM}, a) = 0
f_DE(c::Union{FlatWCDM, ClosedWCDM, OpenWCDM}, a) = -3 * (1 + c.w0) + 3 * c.wa * ((a - 1) / log(a - 1.0e-5) - 1)


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
Cosmology.LCDM{Float64}(0.69, 0.0, 0.7099122024007928, 0.29, 8.77975992071536e-5)

julia> c = cosmology(OmegaK=0.1)
Cosmology.LCDM{Float64}(0.69, 0.1, 0.6099122024007929, 0.29, 8.77975992071536e-5)

julia> c = cosmology(w0=-0.9, OmegaK=-0.1)
Cosmology.WCDM{Float64}(0.69, -0.1, 0.8099122024007929, 0.29, 8.77975992071536e-5, -0.9, 0.0)
```
"""
function cosmology(;
        h = 0.69,
        Neff = 3.04,
        OmegaK = 0,
        OmegaM = nothing,
        OmegaR = nothing,
        OmegaC = 0.25,
        OmegaB = 0.05,
        Tcmb = 2.7255,
        w0 = -1,
        wa = 0
    )

    if !(OmegaR === nothing)
        Tcmb = 0.0
        Neff = 0.0
    end

    if !(OmegaM === nothing)
        OmegaC = OmegaM
        OmegaB = 0.0
    end


    if !(w0 == -1 && wa == 0)
        return WCDM(; h, Ω_k = OmegaK, Ω_c = OmegaC, Ω_b = OmegaB, Neff, Tcmb, w0, wa)
    else
        return LCDM(; h, Ω_k = OmegaK, Ω_c = OmegaC, Ω_b = OmegaB, Neff, Tcmb, w0, wa)
    end

end

# Hubble rate

"""
    scale_factor(z)

Return the scale factor ``a(t)`` for a given redshift ``z(t)``. According to the
[Friedmann–Lemaître–Robertson–Walker metric](https://en.wikipedia.org/wiki/Friedmann–Lemaître–Robertson–Walker_metric)
it's given as ``a = 1/(1 + z)`` ([Schneider 2015, p. 186](@cite Schneider2015)).

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
([Schneider 2015, p. 183](@cite Schneider2015)).
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
function Z end
Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) =
    quadgk(a -> 1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) =
    quadgk(a -> 1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]

@doc raw"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Comoving radial distance (``D_C``) in Mpc at redshift `z₂` as seen by an observer at `z₁`.
Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if
provided.

It's calculated as ``D_C = D_{H0} Z``, where ``D_{H0}`` is the Hubble distance at
the present epoch and, ``Z = \int_{z_1}^{z_2} \frac{dz}{E(z)}``.
"""
function comoving_radial_dist end
comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_transverse_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Comoving transverse distance (``D_C``) in Mpc at redshift `z₂` as seen by an observer at `z₁`.
Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if
provided.

It's identical to the comoving radial distance for a flat cosmological model.

### See also
[`comoving_radial_dist`](@ref)
"""
function comoving_transverse_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...)
    if c.Ω_k > 0
        sqrtΩk = sqrt(c.Ω_k)
        return hubble_dist0(c) * sinh(sqrtΩk * Z(c, z₁, z₂; kws...)) / sqrtΩk
    elseif c.Ω_k < 0
        sqrtΩk = sqrt(abs(c.Ω_k))
        return hubble_dist0(c) * sin(sqrtΩk * Z(c, z₁, z₂; kws...)) / sqrtΩk
    else
        return comoving_radial_dist(c, z₁, z₂; kws...)
    end
end

"""
    angular_diameter_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Ratio of the proper transverse size in Mpc of an object at redshift `z₂` to its
angular size in radians, as seen by an observer at `z₁`.  Redshift `z₁` defaults
to 0 if omitted.  Will convert to compatible unit `u` if provided.
"""
function angular_diameter_dist end
angular_diameter_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) / (1 + z)
angular_diameter_dist(c::AbstractCosmology, z₁, z₂; kws...) =
    comoving_transverse_dist(c, z₁, z₂; kws...) / (1 + z₂)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to
compatible unit `u` if provided.
"""
function luminosity_dist end
luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    distmod(c::AbstractCosmology, z)

Distance modulus in magnitudes at redshift `z`.
"""
function distmod end
distmod(c::AbstractCosmology, z; kws...) =
    5 * log10(luminosity_dist(c, z; kws...) / Mpc) + 25

# volumes

"""
    comoving_volume([u::Unitlike,] c::AbstractCosmology, z)

Comoving volume in cubic Gpc out to redshift `z`. Will convert to compatible unit `u` if provided.
"""
function comoving_volume(c::AbstractCosmology, z; kws...)
    if c.Ω_k == 0
        return (4pi / 3) * (comoving_radial_dist(Gpc, c, z; kws...))^3
    end
    DH = hubble_dist0(Gpc, c)
    x = comoving_transverse_dist(Gpc, c, z; kws...) / DH
    if c.Ω_k > 0
        sqrtΩk = sqrt(c.Ω_k)
        return 2pi * DH^3 * (x * sqrt(1 + c.Ω_k * x^2) - asinh(sqrtΩk * x) / sqrtΩk) / c.Ω_k
    else # c.Ω_k < 0
        sqrtΩk = sqrt(abs(c.Ω_k))
        return 2pi * DH^3 * (x * sqrt(1 + c.Ω_k * x^2) - asin(sqrtΩk * x) / sqrtΩk) / c.Ω_k
    end
end

"""
    comoving_volume_element([u::Unitlike,] c::AbstractCosmology, z)

Comoving volume element in Gpc out to redshift `z`. Will convert to compatible unit `u` if provided.
"""
function comoving_volume_element end
comoving_volume_element(c::AbstractCosmology, z; kws...) =
    hubble_dist0(Gpc, c) * angular_diameter_dist(Gpc, c, z; kws...)^2 / a2E(c, scale_factor(z))

# times

@doc raw"""
    T(c::AbstractCosmology, a0, a1; kws...)

The line-of-sight contributions for lookback time.

It performs the integral
```math
T = \int_{a_0}^{a_1} \frac{da}{a E(a)}
```
"""
T(c::AbstractCosmology, a0, a1; kws...) = quadgk(a -> a / a2E(c, a), a0, a1; kws...)[1]

"""
    age([u::Unitlike,] c::AbstractCosmology, z)

Age of the universe in Gyr at redshift `z`. Will convert to compatible unit `u` if provided.
"""
function age end
age(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, 0, scale_factor(z); kws...)

"""
    lookback_time([u::Unitlike,] c::AbstractCosmology, z)

Difference between age at redshift 0 and age at redshift `z` in Gyr.
Will convert to compatible unit `u` if provided.
"""
function lookback_time end
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, scale_factor(z), 1; kws...)

# Easily select a different unit
for f in (
        :hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time,
        :comoving_radial_dist, :comoving_transverse_dist,
        :angular_diameter_dist, :luminosity_dist,
        :comoving_volume, :comoving_volume_element,
        :age, :lookback_time,
    )
    @eval $f(u::Unitful.Unitlike, args...; kws...) = uconvert(u, $f(args...; kws...))
end

# Density evolution

OmegaM(c::AbstractCosmology, z) = c.Ω_m * (1 + z)^3 / E(c, z)^2
OmegaDE(c::AbstractCosmology, z) = c.Ω_Λ * scale_factor(z)^f_DE(c, scale_factor(z)) / E(c, z)^2


# Growth of structure

function growth_derivatives!(du, u, c, a)
    # Ex 1.118 in Leclerq's thesis shows the equation
    # satisfied by the second order growth factor

    D_1, D_1′, D_2, D_2′ = u
    z = 1 / a - 1
    Omega_m = OmegaM(c, z)
    Omega_de = OmegaDE(c, z)

    D_1′′ = 1.5 * Omega_m * D_1 / a^2 - (Omega_de - 0.5 * Omega_m + 2) * D_1′ / a
    D_2′′ = 1.5 * Omega_m * (D_2 - D_1^2) / a^2 - (Omega_de - 0.5 * Omega_m + 2) * D_2′ / a

    du[1] = D_1′
    du[2] = D_1′′
    du[3] = D_2′
    du[4] = D_2′′
    return du
end

"""
Returns growth factors and rates at scale factor a
Returns: D1(a), D1′(a), D1''(a), D2(a), D2′(a), D2''(a)
"""
function growth_factor(c::AbstractCosmology, a::Vector{<:Real})
    a0 = 1.0e-2
    aspan = (a0, 1.0)
    a_save = a
    push!(a_save, 1.0)
    D1_in = a0 # EdS conditions
    D1′_in = 1.0
    D2_in = -3.0 * D1_in^2 / 7 # Approx eq. 1.119 at tau = 0, Eds implies Om = 1
    D2′_in = -6 * D1_in / 7
    u0 = [D1_in, D1′_in, D2_in, D2′_in]
    prob = ODEProblem(growth_derivatives!, u0, aspan, c)
    sol = solve(prob, saveat = a_save)
    D1 = sol[1, :] ./ sol[1, end]
    D1′ = sol[2, :] ./ sol[1, end]
    D2 = sol[3, :] ./ sol[3, end]
    D2′ = sol[4, :] ./ sol[1, end]
    du = similar(u0)
    D1′′ = similar(a_save)
    D2′′ = similar(a_save)
    for i in eachindex(a_save)
        growth_derivatives!(du, sol[:, i], c, a_save[i])
        D1′′[i] = du[2] / sol[1, end]
        D2′′[i] = du[4] / sol[3, end]
    end
    return D1, D1′, D1′′, D2, D2′, D2′′
end
end # module
