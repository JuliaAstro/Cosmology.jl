module UnitfulExt

using Cosmology
using Unitful, UnitfulAstro
import Unitful: km, s, Gyr
using UnitfulAstro: Mpc, Gpc

# Easily select a different unit
for (f, base_units) in (
        (:hubble_dist0, u"Mpc"),
        (:hubble_dist, u"Mpc"),
        (:hubble_time0, u"Gyr" ),
        (:hubble_time, u"Gyr" ),
        (:comoving_radial_dist, u"Mpc"),
        (:comoving_transverse_dist, u"Mpc"),
        (:angular_diameter_dist, u"Mpc"),
        (:luminosity_dist, u"Mpc"),
        (:comoving_volume, u"Mpc^3"),
        (:age, u"Gyr"),
        (:lookback_time, u"Gyr"),
    )
    @eval Cosmology.$f(u::Unitful.Unitlike, args...; kws...) = uconvert(u, $f(args...; kws...) * $base_units)
end

end # module
