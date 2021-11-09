__precompile__()
module GadgetReader

using LinearAlgebra: rotate!

using Cosmology
using GadgetIO
using LinearAlgebra
using Statistics
using Unitful
using UnitfulAstro

include("snapshot.jl")
include("particles.jl")
include("read_snapshot.jl")
include("transformations.jl")
include("units.jl")
include("utils.jl")

# structs
export Snapshot,
       Particles,

       particle_type_id,

       # read
       read_redshift,
       read_header_particle_mass,

       # transformations
       rotate,
       rotate!,
       rotate_edgeon,
       rotate_edgeon!,
       translate,
       translate!,
       center_of_mass_iterative,
       translate_to_center_of_mass_iterative,
       translate_to_center_of_mass_iterative!,
       center_of_velocity,
       translate_to_center_of_velocity,
       translate_to_center_of_velocity!,


       # units

end
