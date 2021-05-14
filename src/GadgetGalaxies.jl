__precompile__()
module GadgetGalaxies

using LinearAlgebra: rotate!

using Cosmology
using GadgetIO
using LinearAlgebra
using Statistics
using Unitful
using UnitfulAstro

include("snapshot.jl")
include("galaxy.jl")
include("read_halo.jl")
include("transformations.jl")
include("shapes.jl")
include("units.jl")
include("utils.jl")

# structs
export Snapshot,
       Galaxy,
       GalaxyGroup,

       particle_type_id,

       # read halo
       read_halo!,
       read_redshift,
       read_header_particle_mass,
       read_galaxy_prop,
       read_galaxy_pos,
       read_galaxy_vel,
       is_main_halo,

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

       # shapes
       Ellipsoid,
       Sphere,
       rotation_matrix_edgeon,
       rotation_matrix_2D,
       rotation_matrix_axis_ratios,
       r²_ellipsoid,
       triaxiality,
       ellipticity,
       eccentricity,

       # units
       convert_units!,
       convert_units_full,
       convert_units_physical,
       convert_units_physical!,

       convert_units_pos,
       convert_units_full_pos,
       convert_units_physical_pos,
       convert_units_physical_pos!,
       convert_units_vel,
       convert_units_full_vel,
       convert_units_physical_vel,
       convert_units_physical_vel!,
       convert_units_temp,
       convert_units_full_temp,
       convert_units_physical_temp,
       convert_units_physical_temp!,
       convert_units_mass,
       convert_units_full_mass,
       convert_units_physical_mass,
       convert_units_physical_mass!,
       convert_units_age,
       convert_units_full_age,
       convert_units_physical_age,
       convert_units_physical_age!,

       convert_units_solar_metallicity

end
