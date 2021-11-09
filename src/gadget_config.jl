"""
    GadgetConfig(snapprops, subprops, particle_types)

Returns a GADGET Configuration struct with the snap and sub file properties and the mapping
from particle type to the index in the GADGET files.
"""
struct GadgetConfig
    snapprops::Dict{String}
    subprops::Dict{String}
    particle_types::Dict{Symbol,Int}
end

const _snapprops_default, _subprops_default, _particle_types_default = let
    mass = (h -> 1e10 / h.h0, u"Msun")
    pos = (h -> 1 / (h.h0 * (h.z + 1)), u"kpc")
    vel = (h -> 1 / sqrt(h.z + 1), u"km/s")

    snapprops = Dict("POS" => pos, "VEL" => vel, "MASS" => mass)

    subprops = Dict(
        "MTOT" => mass,
        "GPOS" => pos,
        "MVIR" => mass,
        "RVIR" => pos,
        "MSUB" => mass,
        "SPOS" => pos,
        "SVEL" => (h -> 1, u"km/s"), # physical units
    )

    particle_types = Dict(:gas => 0, :dm => 1, :stars => 4, :bh => 5)

    return snapprops, subprops, particle_types
end

"""
    GadgetConfig(; snapprops, subprops, particle_types)

Convenience method for creating a `GadgetConfig` where the properties of the snap or sub file are
not specified.

The default values are as follows:

```julia
using Unitful, UnitfulAstro

mass = (h -> 1e10 / h.h0, u"Msun")
pos = (h -> 1 / (h.h0 * (h.z + 1)), u"kpc")
vel = (h -> 1 / sqrt(h.z + 1), u"km/s")

snapprops = Dict("POS" => pos, "VEL" => vel, "MASS" => mass)

subprops = Dict(
    "MTOT" => mass,
    "GPOS" => pos,
    "MVIR" => mass,
    "RVIR" => pos,
    "MSUB" => mass,
    "SPOS" => pos,
    "SVEL" => (h -> 1, u"km/s"), # physical units
)

particle_types = Dict(:gas => 0, :dm => 1, :stars => 4, :bh => 5)
```
"""
function GadgetConfig(;
    snapprops=_snapprops_default,
    subprops=_subprops_default,
    particle_types=_particle_types_default,
)
    GadgetConfig(snapprops, subprops, particle_types)
end
