Cosmology calculator for Julia
==============================

[![Build Status](https://img.shields.io/travis/JuliaAstro/Cosmology.jl.svg?style=flat-square&label=build)](https://travis-ci.org/JuliaAstro/Cosmology.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaAstro/Cosmology.jl.svg?style=flat-square)](https://coveralls.io/r/JuliaAstro/Cosmology.jl?branch=master)


Installation
------------

To install the package:

```julia
pkg> add Cosmology
```

Then, to load into your session:

```julia
julia> using Cosmology
```

`Cosmology.jl` uses [`Unitful.jl`](https://github.com/ajkeller34/Unitful.jl) and
[`UnitfulAstro.jl`](https://github.com/JuliaAstro/UnitfulAstro.jl) to handle units.

Cosmological Models
-------------------

First, pick a cosmological model using the `cosmology` function,
which takes the following options:

<table>
  <tr>
    <td>h = 0.69</td>
    <td>Dimensionless Hubble constant</td>
  </tr>
  <tr>
    <td>OmegaK = 0</td>
    <td>Curvature density, Ω<sub>k</sub></td>
  </tr>
  <tr>
    <td>OmegaM = 0.29</td>
    <td>Matter density, Ω<sub>m</sub></td>
  </tr>
  <tr>
    <td>OmegaR = Ω<sub>γ</sub> + Ω<sub>ν</sub></td>
    <td>Radiation density, Ω<sub>r</sub></td>
  </tr>
  <tr>
    <td>Tcmb = 2.7255</td>
    <td>CMB temperature (K), used to compute Ω<sub>γ</sub></td>
  </tr>
  <tr>
    <td>Neff = 3.04</td>
    <td>Effective number of massless neutrino species, used to compute Ω<sub>ν</sub></td>
  </tr>
  <tr>
    <td>w0 = -1</td>
    <td>CPL dark energy equation of state, w = w0 + wa*(1-a)</td>
  </tr>
  <tr>
    <td>wa = 0</td>
    <td>CPL dark energy equation of state, w = w0 + wa*(1-a)</td>
  </tr>
</table>

```julia
julia> using Cosmology

julia> c = cosmology()
FlatLCDM(0.69,0.7099122024007928,0.29,8.779759920715362e-5)

julia> c = cosmology(OmegaK=0.1)
OpenLCDM(0.69,0.1,0.6099122024007929,0.29,8.779759920715362e-5)

julia> c = cosmology(w0=-0.9, OmegaK=-0.1)
ClosedWCDM(0.69,-0.1,0.8099122024007929,0.29,8.779759920715362e-5,-0.9,0.0)
```

Distances
---------

<table>
  <tr>
    <td>angular_diameter_dist(cosmo,&nbsp;z)</td>
    <td>Ratio of an object's proper transverse size (in Mpc) to its angular size (in radians)</td>
  </tr>
  <tr>
    <td>comoving_radial_dist(cosmo,&nbsp;z)</td>
    <td>Comoving radial distance to redshift z, in Mpc</td>
  </tr>
  <tr>
    <td>luminosity_dist(cosmo, z)</td>
    <td>Bolometric luminosity distance, in Mpc</td>
  </tr>
  <tr>
    <td>distmod(cosmo, z)</td>
    <td>Distance modulus, in units of magnitude</td>
  </tr>
</table>

```julia
julia> using Cosmology

julia> c = cosmology(OmegaM=0.26)
FlatLCDM(0.69,0.7399122024007928,0.26,8.779759920715362e-5)

julia> angular_diameter_dist(c, 1.2)
1784.0089227105113 Mpc
```

For each function returning a unitful number, you can specify a different unit
for the result as first argument to the function:

```julia
julia> luminosity_dist(c, 1.5)
11420.338287150073 Mpc

julia> luminosity_dist(u"Gpc", c, 1.5)
11.420338287150074 Gpc
```

Volumes
-------

<table>
  <tr>
    <td>comoving_volume_element(cosmo,&nbsp;z)</td>
    <td>Comoving volume element out to redshift z, in Gpc<sup>3</sup></td>
  </tr>
  <tr>
    <td>comoving_volume(cosmo,&nbsp;z)</td>
    <td>Comoving volume out to redshift z, in Gpc<sup>3</sup></td>
  </tr>
</table>

```julia
julia> using Cosmology

julia> c = cosmology(OmegaM=0.26)
FlatLCDM(0.69,0.7399122024007928,0.26,8.779759920715362e-5)

julia> comoving_volume_element(c, 2.1)
46.74459228888612 Gpc^3

julia> comoving_volume(c, 0.6)
49.3633436631307 Gpc^3

julia> comoving_volume(u"ly^3", c, 0.6)
1.7127035381752994e30 ly^3
```

Times
-----

<table>
  <tr>
    <td>age(cosmo, z)</td>
    <td>Age of the universe at redshift z, in Gyr</td>
  </tr>
  <tr>
    <td>lookback_time(cosmo, z)</td>
    <td>Difference between age at redshift 0 and age at redshift z, in Gyr</td>
  </tr>
</table>

```julia
julia> using Cosmology

julia> c = cosmology(OmegaM=0.26)
FlatLCDM(0.69,0.7399122024007928,0.26,8.779759920715362e-5)

julia> age(c, 1.2)
5.445600787626434 Gyr

julia> lookback_time(u"yr", c, 1.2)
8.761660748088268e9 yr
```
