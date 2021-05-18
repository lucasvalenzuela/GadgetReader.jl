```@meta
CurrentModule = GadgetGalaxies
```

# GadgetGalaxies

[`GadgetGalaxies.jl`](https://github.com/lucasvalenzuela/GadgetGalaxies.jl) provides a wrapper for [`GadgetIO.jl`](https://github.com/LudwigBoess/GadgetIO.jl) to facilitate reading groups and subhalos from subfind and snapshot particle data, adding conversion functionality to physical units and various galactic quantities.


## Installation

The latest version of the package is available for Julia 1.6 and newer versions. It is recommended to install `GadgetIO` alongside it:

```julia
using Pkg
Pkg.add("GadgetIO")
Pkg.add("GadgetGalaxies")
```

## Getting Started

First, we need to create a [`Snapshot`](@ref) (using the filebase names without ".0", ".1", etc.):

```julia
snapshot = Snapshot("filebase/of/snapfiles", "filebase"/of/subfiles")
```

If your simulations use the format `"groups_XXX/sub_XXX"` and `"snapdir_XXX/snap_XXX"`, you can also initialize the snapshot like this, for example for snapshot 140:

```julia
simdir = "directory/of/simulation"
snapshot = Snapshot(simdir, 140)
```

To prepare a [`Galaxy`](@ref) or [`GalaxyGroup`] for reading its particles, pass it the zero-based index of the subhalo in subfind:

```julia
g = Galaxy(snapshot, 1234)
gr = GalaxyGroup(snapshot, 123)
```

After the galaxy has been set up, use [`read_halo!`](@ref) to read in some basic properties of all stellar, dark matter, and gas particles (by default within the half-mass radius), which can be accessed in multiple ways:

```julia
read_halo!(g)
g.gas["POS"]
g.gas.pos # alternative notation
g[:gas].pos # alternative notation
```
