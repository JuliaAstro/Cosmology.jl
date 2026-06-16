```@meta
DocTestSetup = quote
    using Cosmology, DynamicQuantities
end
```

# API/Reference


!!! tip "Units"
    [DynamicQuantities.jl](https://github.com/JuliaPhysics/DynamicQuantities.jl) and [Unitful.jl](https://github.com/JuliaPhysics/Unitful.jl) work seamlessly with Cosmology.jl. In order to use their features, make sure to run either:

    ```julia-repl
    pkg> add DynamicQuantities
    julia> using DynamicQuantities
    ```

    or:

    ```julia-repl
    pkg> add Unitful UnitfulAstro
    julia> using Unitful, UnitfulAstro
    ```


## Cosmological Models

```@docs
cosmology
```

## Distances

```@docs
angular_diameter_dist
comoving_radial_dist
luminosity_dist
distmod
hubble_dist
```

### Examples

```jldoctest
julia> c = cosmology(OmegaM=0.26)
Cosmology.FlatLCDM{Float64}(0.69, 0.7399122024007928, 0.26, 8.77975992071536e-5)

julia> angular_diameter_dist(c, 1.2)
1784.0089227105113 Mpc

julia> angular_diameter_dist(c, 0.7, 1.2)
606.6521737365097 Mpc

julia> luminosity_dist(c, 1.5)
11420.338287150073 Mpc

julia> luminosity_dist(us"Constants.Gpc", c, 1.5) # Can convert to appropriate unit
11.420338287150072 Gpc
```


## Volumes

```@docs
comoving_volume_element
comoving_volume
```

### Examples

```jldoctest
julia> c = cosmology(OmegaM=0.26)
Cosmology.FlatLCDM{Float64}(0.69, 0.7399122024007928, 0.26, 8.77975992071536e-5)

julia> comoving_volume_element(c, 2.1)
46.74459228888613 Gpc³

julia> comoving_volume(c, 0.6)
49.363343663130685 Gpc³

julia> comoving_volume(us"Constants.ly^3", c, 0.6)
1.7127035381752996e30 ly³
```

## Times

```@docs
age
lookback_time
hubble_time
```

### Examples
```jldoctest
julia> c = cosmology(OmegaM=0.26)
Cosmology.FlatLCDM{Float64}(0.69, 0.7399122024007928, 0.26, 8.77975992071536e-5)

julia> age(c, 1.2)
5.4454795007229455 Gyr

julia> lookback_time(us"yr", c, 1.2)
8.761465604385489e9 yr
```

## Miscellaneous
```@docs
H
scale_factor
```


## Bibliography
```@bibliography
Pages = ["api.md"]
```
