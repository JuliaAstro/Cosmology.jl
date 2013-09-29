module Cosmology

export cosmology,
       age_gyr,
       angular_diameter_dist_mpc,
       comoving_radial_dist_mpc,
       comoving_transverse_dist_mpc,
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

immutable FlatLCDM <: AbstractFlatCosmology
    h::Float64
    Ω_Λ::Float64
    Ω_m::Float64
    Ω_r::Float64

    function FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real)
        new(float64(h), float64(Ω_Λ), float64(Ω_m), float64(Ω_r))
    end
end

a2E(c::FlatLCDM, a::Float64) = sqrt(c.Ω_r + c.Ω_m*a + c.Ω_Λ*a^4)

immutable ClosedLCDM <: AbstractClosedCosmology
    h::Float64
    Ω_k::Float64
    Ω_Λ::Float64
    Ω_m::Float64
    Ω_r::Float64

    function ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real)
        new(float64(h), float64(Ω_k), float64(Ω_Λ), float64(Ω_m), float64(Ω_r))
    end
end

immutable OpenLCDM <: AbstractOpenCosmology
    h::Float64
    Ω_k::Float64
    Ω_Λ::Float64
    Ω_m::Float64
    Ω_r::Float64

    function OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real)
        new(float64(h), float64(Ω_k), float64(Ω_Λ), float64(Ω_m), float64(Ω_r))
    end
end

function a2E(c::Union(ClosedLCDM,OpenLCDM), a::Float64)
    a2 = a*a
    sqrt(c.Ω_r + c.Ω_m*a + (c.Ω_k + c.Ω_Λ*a2)*a2)
end

for c in ("Flat", "Open", "Closed")
    @eval begin
        immutable $(symbol("$(c)WCDM")) <: $(symbol("Abstract$(c)Cosmology"))
            h::Float64
            Ω_k::Float64
            Ω_Λ::Float64
            Ω_m::Float64
            Ω_r::Float64
            w0::Float64
            wa::Float64
        end
    end
end

function WCDM(h::Float64, Ω_k::Float64, Ω_Λ::Float64, Ω_m::Float64, Ω_r::Float64, w0::Float64, wa::Float64)
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    end
end

function a2E(c::Union(FlatWCDM,ClosedWCDM,OpenWCDM), a::Float64)
    ade = exp((1.0 - 3.0*(c.w0 + c.wa))*log(a) + 3.0*c.wa*(a - 1.0))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k*a)*a + c.Ω_Λ*ade)
end

function cosmology(;h=0.7,
                   Neff=3.04,
                   OmegaK=0,
                   OmegaM=0.3,
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
        return WCDM(float64(h), float64(OmegaK), float64(OmegaL),
                    float64(OmegaM), float64(OmegaR), float64(w0), float64(wa))
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

# times

T(c::AbstractCosmology, a0::Float64, a1::Float64) = ((q,_) = quadgk(x::Float64->x/a2E(c,x), a0, a1); q)
age_gyr(c::AbstractCosmology, z) = hubble_time_gyr0(c)*T(c, 0., scale_factor(z))
lookback_time_gyr(c::AbstractCosmology, z) = hubble_time_gyr0(c)*T(c, scale_factor(z), 1.)

end # module
