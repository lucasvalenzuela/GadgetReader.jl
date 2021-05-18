```@meta
CurrentModule = GadgetGalaxies
```

# Kinematic Parameters

`GadgetGalaxies` provides various functions to compute kinematic parameters of galaxies.

## Angular Momentum

The (specific) angular momentum can easily be computed from a [`Galaxy`](@ref) or [`GalaxyGroup`](@ref):

```@docs
angular_momentum
```

```@docs
specific_angular_momentum
```

## B-Value

The b-value methods for a galaxy are called in the same way as the angular momenta methods. However, since this is generally only computed for the stellar component, it is recommended to center the stars positionally and kinematically:

```julia
using Unitful, UnitfulAstro

gr::GalaxyGroup # assume an already read-in galaxy group with :full units

rvir = read_galaxy_prop(gr, "RVIR")

# first half-mass radius estimate
r = half_mass_radius(gr.stars; rmax=0.1rvir)

# shift to center
translate_to_center_of_mass_iterative!(gr, 4r, :stars)

# second determination of the half-mass radius (this time with proper centering)
r = half_mass_radius(gr.stars; rmax=0.1rvir)

# center velocity
translate_to_center_of_velocity!(gr, 4r)

# finally determine b_value within two half-mass radii
b = b_value(gr.stars, Sphere(2r))
```

```@docs
b_value
```
