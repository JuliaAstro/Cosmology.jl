module UnitfulExt

using Cosmology
import Unitful as U
import DynamicQuantities as DQ

# Easily select a different unit
for f in (
    :age,
    :angular_diameter_dist,
    :comoving_radial_dist,
    :comoving_transverse_dist,
    :comoving_volume,
    :comoving_volume_element,
    :hubble_dist,
    :hubble_time,
    :luminosity_dist,
    :lookback_time,
    )
    @eval function Cosmology.$f(u::U.Unitlike, args...; kws...)
        q_dq = $f(args...; kws...) |> DQ.uexpand
        q_unitful = convert(U.Quantity, q_dq)
        return U.uconvert(u, q_unitful)
    end
end

end # module
