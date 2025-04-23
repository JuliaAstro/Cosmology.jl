```@meta
CurrentModule = Cosmology
```
# Internals

The following types and methods are internal, and should not be considered safe for public use.

## Types
```@docs
AbstractCosmology
FlatLCDM
ClosedLCDM
OpenLCDM
```

## Methods
```@docs
E
Z
a2E
a2E(::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
hubble_dist0
hubble_time0
```

## Bibliography

```@bibliography
Pages = ["internals.md"]
```
