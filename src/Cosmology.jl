module Cosmology

using QuadGK, Unitful
import Unitful: km, s, Gyr
using UnitfulAstro: Mpc, Gpc
using DifferentialEquations


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
       growth_factor,
       σ_T, m_H, m_He, m_p



abstract type AbstractCosmology end
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end


for model in ("LCDM", "WCDM")
    for curv in ("Flat", "Open", "Closed")
        name = Symbol("$(curv)$(model)")
        @eval begin
            Base.@kwdef struct $(name){N <: Real} <: $(Symbol("Abstract$(curv)Cosmology"))
                h::N = 0.67
                Ω_k::N = 0
                Ω_c::N = 0.3
                Ω_b::N = 0
                Neff::N = 3.04
                T_cmb::N = 2.755
                w0::N = -1
                wa::N = 0

                # Derived densities
                Ω_γ::N = 4.48131e-7 * T_cmb^4 / h^2
                Ω_ν::N = Neff * Ω_γ * (7 / 8) * (4 / 11)^(4 / 3)
                Ω_r::N = Ω_γ + Ω_ν
                Ω_m::N = Ω_c + Ω_b
                Ω_Λ::N = 1. - Ω_m - Ω_r - Ω_k  
                

            end
            function $(name)(args...)
                $(name)(promote(map(float, args)...)...)
            end
            function $(name)(;kwargs...)
                $(name)(;promote(map(float, kwargs)...)...)
            end
            
        end
    end
    @eval begin
        function $(Symbol("$model"))(;kwargs...)
            kwargs = Dict(kwargs)
            Ω_k = haskey(kwargs, :Ω_k) ? kwargs[:Ω_k] : 0.
            if Ω_k < 0
                $(Symbol("Closed$(model)"))(;kwargs...)
            elseif Ω_k > 0
                $(Symbol("Open$(model)"))(;kwargs...)
            else
                $(Symbol("Flat$(model)"))(;kwargs...)
            end
        end
    end
    @eval begin
        function $(Symbol("$model"))(args...)
            Ω_k = args[2]
            if Ω_k < 0
                $(Symbol("Closed$(model)"))(args...)
            elseif Ω_k > 0
                $(Symbol("Open$(model)"))(args...)
            else
                $(Symbol("Flat$(model)"))(args...)
            end
        end
    end

end

function a2E(c::Union{FlatLCDM,ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end
f_DE(c::Union{FlatLCDM,ClosedLCDM,OpenLCDM}, a) = 0

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end
f_DE(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a) = (-3 * (1 + c.w0) + 3 * c.wa * ((a - 1) / log(a - 1e-5) - 1))

"""
    cosmology(;h = 0.69,
               Neff = 3.04,
               OmegaK = 0,
               OmegaM = 0.29,
               OmegaR = nothing,
               Tcmb = 2.7255,
               w0 = -1,
               wa = 0)


# Parameters
* `h` - Dimensionless Hubble constant
* `OmegaK` - Curvature density (Ω_k)
* `OmegaM` - Matter density (Ω_m)
* `OmegaR` - Radiation density (Ω_r)
* `Tcmb` - CMB temperature in Kelvin; used to compute Ω_γ
* `Neff` - Effective number of massless neutrino species; used to compute Ω_ν
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
                   OmegaC = 0.25,
                   OmegaB = 0.05,
                   OmegaR = nothing,
                   OmegaM = nothing,
                   Tcmb = 2.7255,
                   w0 = -1,
                   wa = 0)
    
    if !(OmegaR === nothing)
        Tcmb = 0.
        Neff = 0.
    end

    if !(OmegaM === nothing)
        OmegaC = OmegaM
        OmegaB = 0.
    end


    if !(w0 == -1 && wa == 0)
        return WCDM(;h = h, Ω_k = OmegaK, Ω_c = OmegaC, Ω_b = OmegaB, Neff = Neff, T_cmb = Tcmb, w0 = w0, wa = wa)
    else
        return LCDM(;h = h, Ω_k = OmegaK, Ω_c = OmegaC, Ω_b = OmegaB, Neff = Neff, T_cmb = Tcmb, w0 = w0, wa = wa)
    end

    
end

# hubble rate

scale_factor(z) = 1 / (1 + z)
E(c::AbstractCosmology, z) = (a = scale_factor(z); a2E(c, a) / a^2)
H(c::AbstractCosmology, z) = 100 * c.h * E(c, z) * km / s / Mpc

hubble_dist0(c::AbstractCosmology) = 2997.92458 / c.h * Mpc
hubble_dist(c::AbstractCosmology, z) = hubble_dist0(c) / E(c, z)

hubble_time0(c::AbstractCosmology) = 9.777922216807891 / c.h * Gyr
hubble_time(c::AbstractCosmology, z) = hubble_time0(c) / E(c, z)

# distances

Z(c::AbstractCosmology, z::Real, ::Nothing; kws...) =
    QuadGK.quadgk(a->1 / a2E(c, a), scale_factor(z), 1; kws...)[1]
Z(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) =
    QuadGK.quadgk(a->1 / a2E(c, a), scale_factor(z₂), scale_factor(z₁); kws...)[1]

comoving_radial_dist(c::AbstractCosmology, z₁, z₂ = nothing; kws...) = hubble_dist0(c) * Z(c, z₁, z₂; kws...)

"""
    comoving_radial_dist([u::Unitlike,] c::AbstractCosmology, [z₁,] z₂)

Comoving radial distance in Mpc at redshift `z₂` as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided.
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

Ratio of the proper transverse size in Mpc of an object at redshift `z₂` to its angular size in radians, as seen by an observer at `z₁`.  Redshift `z₁` defaults to 0 if omitted.  Will convert to compatible unit `u` if provided.
"""
angular_diameter_dist

luminosity_dist(c::AbstractCosmology, z; kws...) =
    comoving_transverse_dist(c, z; kws...) * (1 + z)

"""
    luminosity_dist([u::Unitlike,] c::AbstractCosmology, z)

Bolometric luminosity distance in Mpc at redshift `z`. Will convert to compatible unit `u` if provided.
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

T(c::AbstractCosmology, a0, a1; kws...) = QuadGK.quadgk(x->x / a2E(c, x), a0, a1; kws...)[1]
"""
    age([u::Unitlike,] c::AbstractCosmology, z)

Age of the universe in Gyr at redshift `z`. Will convert to compatible unit `u` if provided.
"""
age(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, 0, scale_factor(z); kws...)

"""
    lookback_time([u::Unitlike,] c::AbstractCosmology, z)

Difference between age at redshift 0 and age at redshift `z` in Gyr. Will convert to compatible unit `u` if provided.
"""
lookback_time(c::AbstractCosmology, z; kws...) = hubble_time0(c) * T(c, scale_factor(z), 1; kws...)

# Easily select a different unit
for f in (:hubble_dist0, :hubble_dist, :hubble_time0, :hubble_time, :comoving_radial_dist,
          :comoving_transverse_dist, :angular_diameter_dist, :luminosity_dist,
          :comoving_volume, :comoving_volume_element, :age, :lookback_time)
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

function growth_factor(c::AbstractCosmology, a::Vector{<:Real})
    """
    Returns growth factors and rates at scale factor a
    Returns: D1(a), D1′(a), D1''(a), D2(a), D2′(a), D2''(a)
    """
    a0 = 1e-2
    aspan = (a0,1.)
    a_save = a
    push!(a_save, 1.)
    D1_in = a0 # EdS conditions
    D1′_in = 1.
    D2_in = -3. * D1_in^2 / 7 # Approx eq. 1.119 at tau = 0, Eds implies Om = 1
    D2′_in = -6* D1_in / 7
    u0 = [D1_in, D1′_in, D2_in, D2′_in]
    prob = ODEProblem(growth_derivatives!,u0,aspan,c)
    sol = solve(prob, saveat = a_save)
    D1 = sol[1,:] ./ sol[1,end]
    D1′ = sol[2,:] ./ sol[1,end]
    D2 = sol[3,:] ./ sol[3,end]
    D2′ = sol[4,:] ./ sol[1,end]
    du = similar(u0)
    D1′′ = similar(a_save)
    D2′′ = similar(a_save)
    for i in eachindex(a_save)
        growth_derivatives!(du, sol[:,i], c, a_save[i])
        D1′′[i] = du[2] / sol[1,end]
        D2′′[i] = du[4] / sol[3,end]
    end
    return D1, D1′, D1′′, D2, D2′, D2′′
end



end # module

