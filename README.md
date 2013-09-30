Cosmology calculator for Julia
==============================

Installation
------------

To install the package:

    julia> Pkg.add("Cosmology")

Then, to load into your session:

    julia> using Cosmology

Cosmological Models
-------------------

First, pick a cosmological model using the `cosmology` function,
which takes the following options:

<table>
  <tr>
    <td>h = 0.7</td>
    <td>Dimensionless Hubble constant</td>
  </tr>
  <tr>
    <td>OmegaK = 0</td>
    <td>Curvature density, Ω<sub>k</sub></td>
  </tr>
  <tr>
    <td>OmegaM = 0.3</td>
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

```jlcon
julia> c = cosmology()
FlatLCDM(0.7,0.6999146929857499,0.3,8.530701425005273e-5)

julia> c = cosmology(OmegaK=0.1)
OpenLCDM(0.7,0.1,0.5999146929857501,0.3,8.530701425005273e-5)
```

Distances
---------

<table>
  <tr>
    <td>angular_diameter_dist_mpc(cosmo,&nbsp;z)</td>
    <td>Ratio of an object's proper transverse size (in Mpc) to its angular size (in radians)</td>
  </tr>
  <tr>
    <td>comoving_radial_dist_mpc(cosmo,&nbsp;z)</td>
    <td>Comoving radial distance to redshift z, in Mpc</td>
  </tr>
  <tr>
    <td>comoving_volume_gpc3(cosmo,&nbsp;z)</td>
    <td>Comoving volume out to redshift z, in Gpc<sup>3</sup></td>
  </tr>
  <tr>
    <td>luminosity_dist_mpc(cosmo, z)</td>
    <td>Bolometric luminosity distance, in Mpc</td>
  </tr>
</table>

```jlcon
julia> using Cosmology

julia> c = cosmology(OmegaM=0.26)
FlatLCDM(0.7,0.739914695489689,0.26,8.530451031095114e-5)

julia> angular_diameter_dist_mpc(c, 1.2)
1758.5291281199122
```

Times
-----

<table>
  <tr>
    <td>age_gyr(cosmo, z)</td>
    <td>Age of the universe at redshift z, in Gyr</td>
  </tr>
  <tr>
    <td>lookback_time_gyr(cosmo, z)</td>
    <td>Difference between age at redshift 0 and age at redshift z, in Gyr</td>
  </tr>
</table>

```jlcon
julia> using Cosmology

julia> c = cosmology(OmegaM=0.26)
FlatLCDM(0.7,0.739914695489689,0.26,8.530451031095114e-5)

julia> age_gyr(c, 1.2)
5.367964753127867
```

