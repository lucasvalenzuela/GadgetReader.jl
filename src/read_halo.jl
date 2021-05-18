"""
    shift_across_box_border(x::Number, x₀::Number, boxsize::Number, boxsize_half::Number)

Shift coordinate `x` across the box border if reference coordinate `x₀` is on the other side.
"""
function shift_across_box_border(x::Number, x₀::Number, boxsize::Number, boxsize_half::Number)
    if x - x₀ > boxsize_half
        return x - boxsize
    elseif x₀ - x > boxsize_half
        return x + boxsize
    end
    return x
end

"""
    read_halo!(g::Galaxy; kwargs...)

Reads particles of a [`Galaxy`](@ref).

# Keywords
- `radius::Union{Number,Nothing}=nothing`: if `nothing`, the half-mass radius "RHMS"
- `radius_units::Symbol=:physical`: `:physical` or `:sim`, depending on the units of `radius`
- `rad_scale::Real=1`: read particles within `rad_scale * radius` of halo position
- `units::Symbol=:full`: use `:full` for units as `Unitful` quantities, `:physical` for values converted
  to physical units (kpc, km/s, solar metallicities etc.), `:sim` for values in simulation units
- `use_keys::Bool=true`: if Peano-Hilbert keys should be used if possible
- `verbose::Bool=false`: verbose output
- `props::Tuple=(
       (:gas, ["POS", "VEL", "MASS", "TEMP", "SFR", "Zs"]),
       (:dm, ["POS", "VEL"]),
       (:stars, ["POS", "VEL", "MASS", "AGE", "Zs", "iM"])
   )`: particle types (see [`Particles`](@ref)) and properties to be read
"""
function read_halo!(
    g::AbstractGalaxy;
    radius::Union{Number,Nothing}=nothing,
    radius_units::Symbol=:physical,
    rad_scale::Real=1,
    units::Symbol=:full,
    use_keys::Bool=true,
    verbose::Bool=false,
    props::Tuple=(
        (:gas, ["POS", "VEL", "MASS", "TEMP", "SFR", "Zs"]),
        (:dm, ["POS", "VEL"]),
        (:stars, ["POS", "VEL", "MASS", "AGE", "Zs", "iM"]),
    ),
)
    snapbase = string(g.snapshot.snapbase)
    subbase = string(g.snapshot.subbase)

    radius_sim = (isnothing(radius) || radius_units === :sim) ? radius : simulation_units_pos(radius)

    # get global halo properties in simulation units
    halo_pos = read_galaxy_pos(g, :sim; verbose)
    halo_vel = read_galaxy_vel(g, :sim; verbose)

    # Provides properties `h0`, `z`, `time`, `omega_0`, `omega_l`, `num_files`
    h = read_header(snapbase)

    # handle particles for different types
    Threads.@threads for (ptype, blocks) in props
        ptype_id = particle_type_id(ptype)

        prop_data = read_particles_in_halo(
            snapbase,
            blocks,
            subbase,
            getid(g);
            parttype=ptype_id,
            halo_type=get_subfind_type(g),
            # radius=radius_sim,
            rad_scale,
            use_keys,
            verbose,
        )
        p = Particles(ptype, prop_data)

        # shift across box borders and transform into halo frame of reference
        if "POS" in blocks && length(p.pos) > 0
            p.pos .= shift_across_box_border.(p.pos, halo_pos, h.boxsize, 1 // 2 * h.boxsize) .- halo_pos
        end

        # transform velocities into halo frame of reference
        if "VEL" in blocks && length(p.vel) > 0
            p.vel .-= halo_vel
        end

        convert_units!(p, h, units)

        # read particle mass if available in header mass array
        if "MASS" ∉ keys(p)
            p_mass = read_header_particle_mass(g.snapshot, ptype, units)
            # ustrip removes any potential unit to be comparable with 0 (Int64)
            ustrip(p_mass) > 0 && (p.mass = p_mass)
        end

        # add particles to galaxy
        g[ptype] = p
    end

    return g
end


"""
    read_redshift(snapshot::Snapshot)

Returns the redshift ``z`` from the first subfile's header.
"""
function read_redshift(snapshot::Snapshot)
    h = read_subfind_header(snapshot.subbase |> string)
    return h.z
end

"""
    read_dm_particle_mass(snapshot::Snapshot, ptype::Symbol[, units=:full]; verbose::Bool=false)

Returns the particle mass from the first snapfile's header for a particle type (`:dm`, `:gas`, etc.).
"""
function read_header_particle_mass(
    snapshot::Snapshot,
    ptype::Symbol,
    units::Symbol=:full;
    verbose::Bool=false,
)
    h = read_header(snapshot.snapbase |> string)
    dm_mass = h.massarr[particle_type_id(ptype) + 1]
    return convert_units_mass(dm_mass, h, units)
end

"""
    read_galaxy_prop(g::AbstractGalaxy, prop::AbstractString, units::Symbol=:full; verbose::Bool=false)

Reads galaxy property from subfind for `g`. Use `:full` for `units` as `Unitful` quantities,
`:physical` for values converted, and `:sim` for values in simulation units.
"""
function read_galaxy_prop(g::AbstractGalaxy, prop::AbstractString, units::Symbol=:full; verbose::Bool=false)
    subbase = string(g.snapshot.subbase)

    val = read_halo_prop(subbase, getid(g), prop; verbose)

    # if units are not supposed to be converted
    units === :sim && return val

    # get header of snapshot for conversion information
    h = read_subfind_header(subbase)
    return convert_units_subfind_prop(val, prop, h, units; verbose)
end

"""
    read_galaxy_pos(g::AbstractGalaxy, units::Symbol=:full; verbose::Bool=false)

Returns the galaxy's position from subfind using GPOS for a group and SPOS for a subhalo.
"""
function read_galaxy_pos(g::Galaxy, units::Symbol=:full; verbose::Bool=false)
    read_galaxy_prop(g, "SPOS", units; verbose)
end
function read_galaxy_pos(g::GalaxyGroup, units::Symbol=:full; verbose::Bool=false)
    read_galaxy_prop(g, "GPOS", units; verbose)
end


"""
    read_galaxy_vel(g::Galaxy, units::Symbol=:full; verbose::Bool=false)

Returns the galaxy's velocity from subfind using SVEL of FSUB for a group and SVEL for a subhalo.
"""
function read_galaxy_vel(g::Galaxy, units::Symbol=:full; verbose::Bool=false)
    read_galaxy_prop(g, "SVEL", units; verbose)
end
function read_galaxy_vel(g::GalaxyGroup, units::Symbol=:full; verbose::Bool=false)
    subbase = string(g.snapshot.subbase)
    ifsub = read_halo_prop(subbase, g.groupid, "FSUB"; verbose)
    nfiles = read_subfind_header(subbase).num_files
    val, _ = read_halo_prop_and_id(subbase, ifsub, "SVEL", nfiles; verbose)

    # if units are not supposed to be converted
    units === :sim && return val

    h = read_subfind_header(subbase)
    return convert_units_subfind_prop(val, "SVEL", h, units)
end

"""
    is_main_halo(g::Galaxy)

Returns if the subhalo is the first subhalo in its respective group.
"""
function is_main_halo(g::Galaxy)
    subbase = string(g.snapshot.subbase)
    nfiles = read_subfind_header(subbase).num_files
    grnr = read_halo_prop(subbase, g.subid, "GRNR"; verbose=false)
    ifsub, _ = read_halo_prop_and_id(subbase, grnr, "FSUB", nfiles; verbose=false)
    return g.isub == ifsub
end

"""
    get_group(g::Galaxy)

Returns the corresponding [`GalaxyGroup`](@ref) of the given [`Galaxy`](@ref).
"""
function get_group(g::Galaxy)
    subbase = string(g.snapshot.subbase)
    grnr = read_halo_prop(subbase, g.subid, "GRNR"; verbose=false)
    return GalaxyGroup(g.snapshot, grnr)
end

"""
    get_first_subhalo(gr::GalaxyGroup)

Returns the first Subhalo [`Galaxy`](@ref) of the given [`GalaxyGroup`](@ref).
"""
function get_first_subhalo(gr::GalaxyGroup)
    subbase = string(gr.snapshot.subbase)
    ifsub = read_halo_prop(subbase, gr.groupid, "FSUB"; verbose=false)
    return Galaxy(gr.snapshot, ifsub)
end
