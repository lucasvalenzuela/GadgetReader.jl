__precompile__()
module GadgetGalaxies

using Cosmology
using GadgetIO
using Unitful
using UnitfulAstro

include("snapshot.jl")
include("galaxy.jl")
include("read_halo.jl")
include("units.jl")

# structs
export Snapshot,
       Galaxy,
       GalaxyGroup,

       particle_type_id,

       # read halo
       read_halo!,
       read_redshift,

       # units
       convert_units!,
       convert_units_full,
       convert_units_physical,
       convert_units_physical!,

       convert_units_full_pos,
       convert_units_physical_pos,
       convert_units_physical_pos!,
       convert_units_full_vel,
       convert_units_physical_vel,
       convert_units_physical_vel!,
       convert_units_full_temp,
       convert_units_physical_temp,
       convert_units_physical_temp!,
       convert_units_full_mass,
       convert_units_physical_mass,
       convert_units_physical_mass!,
       convert_units_full_age,
       convert_units_physical_age,
       convert_units_physical_age!,

       convert_units_solar_metallicity

end
