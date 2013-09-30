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

```jlcon
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
FlatLCDM(0.69,0.7399122024007928,0.26,8.779759920715362e-5)

julia> angular_diameter_dist_mpc(c, 1.2)
1784.0089227105118
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
FlatLCDM(0.69,0.7399122024007928,0.26,8.779759920715362e-5)

julia> age_gyr(c, 1.2)
5.445600787626434
```

