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
include("global_parameters.jl")
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
       center_of_velocity,
       translate_to_center_of_velocity,
       translate_to_center_of_velocity!,

       # shapes
       Ellipsoid,
       Sphere,
       rotation_matrix_edgeon,
       rotation_matrix_2D,
       rotation_matrix_axis_ratios,
       r²_sphere,
       r²_circle,
       r²_ellipsoid,
       triaxiality,
       ellipticity,
       eccentricity,

       # global parameters
       half_mass_radius,
       half_mass_radius_2D,

       # units
       convert_units!,
       convert_units_full,
       convert_units_physical,
       convert_units_physical!,

       simulation_units_pos,
       simulation_units_pos!,
       convert_units_pos,
       convert_units_full_pos,
       convert_units_physical_pos,
       convert_units_physical_pos!,
       simulation_units_vel,
       simulation_units_vel!,
       convert_units_vel,
       convert_units_full_vel,
       convert_units_physical_vel,
       convert_units_physical_vel!,
       simulation_units_temp,
       simulation_units_temp!,
       convert_units_temp,
       convert_units_full_temp,
       convert_units_physical_temp,
       convert_units_physical_temp!,
       simulation_units_mass,
       simulation_units_mass!,
       convert_units_mass,
       convert_units_full_mass,
       convert_units_physical_mass,
       convert_units_physical_mass!,
       simulation_units_age,
       simulation_units_age!,
       convert_units_age,
       convert_units_full_age,
       convert_units_physical_age,
       convert_units_physical_age!,

       convert_units_solar_metallicity

end
