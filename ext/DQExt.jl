module DQExt

using Cosmology
using DynamicQuantities: Quantity, uconvert
using DynamicQuantities.SymbolicConstants: Mpc
using DynamicQuantities.SymbolicUnits: Gyr

# Easily select a different unit
for (f, base_units) in (
        (:hubble_dist0, Mpc),
        (:hubble_dist, Mpc),
        (:hubble_time0, Gyr ),
        (:hubble_time, Gyr ),
        (:comoving_radial_dist, Mpc),
        (:comoving_transverse_dist, Mpc),
        (:angular_diameter_dist, Mpc),
        (:luminosity_dist, Mpc),
        (:comoving_volume, Mpc^3),
        (:comoving_volume_element, Mpc^3),
        (:age, Gyr),
        (:lookback_time, Gyr),
    )
    @eval Cosmology.$f(u::Quantity, args...; kws...) = uconvert(u, $f(args...; kws...) * $base_units)
end

end # module
