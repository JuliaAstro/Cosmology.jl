```@meta
DocTestSetup = quote
    using Cosmology, Unitful, UnitfulAstro
    ENV["UNITFUL_FANCY_EXPONENTS"] = false
end
```

# API/Reference


!!! tip "Unitful"
    [Unitful.jl](https://github.com/painterqubits/Unitful.jl) works seamlessly with Cosmology.jl. In order to use its features, make sure it is installed and imported, along with [UnitfulAstro](https://github.com/juliaastro/UnitfulAstro.jl).
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

julia> luminosity_dist(u"Gpc", c, 1.5) # Can convert to appropriate unit
11.420338287150074 Gpc
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
46.74459228888613 Gpc^3

julia> comoving_volume(c, 0.6)
49.3633436631307 Gpc^3

julia> comoving_volume(u"ly^3", c, 0.6)
1.7127035381753e30 ly^3
```

## Times

```@docs
age
lookback_time
```

### Examples
```jldoctest
julia> c = cosmology(OmegaM=0.26)
Cosmology.FlatLCDM{Float64}(0.69, 0.7399122024007928, 0.26, 8.77975992071536e-5)

julia> age(c, 1.2)
5.4454795007229455 Gyr

julia> lookback_time(u"yr", c, 1.2)
8.761465604385489e9 yr
```



## Bibliography
```@bibliography
Pages = ["api.md"]
```
