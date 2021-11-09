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
# function read_halo!(
    # g::AbstractGalaxy;
    # radius::Union{Number,Nothing}=nothing,
    # radius_units::Symbol=:physical,
    # rad_scale::Real=1,
    # units::Symbol=:full,
    # use_keys::Bool=true,
    # verbose::Bool=false,
    # props::Tuple=(
        # (:gas, ["POS", "VEL", "MASS", "TEMP", "SFR", "Zs"]),
        # (:dm, ["POS", "VEL"]),
        # (:stars, ["POS", "VEL", "MASS", "AGE", "Zs", "iM"]),
    # ),
# )
    # snapbase = string(g.snapshot.snapbase)
    # subbase = string(g.snapshot.subbase)

    # # Provides properties `h0`, `z`, `time`, `omega_0`, `omega_l`, `num_files`
    # h = read_header(snapbase)

    # # use 1.0radius to convert to float to prevent simulation_units_pos from trying
    # # to convert to Int if radius is int
    # radius_sim = (isnothing(radius) || radius_units === :sim) ? radius : simulation_units_pos(1.0radius, h)

    # # get global halo properties in simulation units
    # halo_pos = read_galaxy_pos(g, :sim; verbose)
    # halo_vel = read_galaxy_vel(g, :sim; verbose)

    # # handle particles for different types
    # Threads.@threads for (ptype, blocks) in props
        # ptype_id = particle_type_id(ptype)

        # prop_data = read_particles_in_halo(
            # snapbase,
            # blocks,
            # subbase,
            # getid(g);
            # parttype=ptype_id,
            # halo_type=get_subfind_type(g),
            # radius=radius_sim,
            # rad_scale,
            # use_keys,
            # verbose,
        # )
        # p = Particles(ptype, prop_data)

        # # shift across box borders and transform into halo frame of reference
        # if "POS" in blocks && length(p.pos) > 0
            # p.pos .= shift_across_box_border.(p.pos, halo_pos, h.boxsize, 1 // 2 * h.boxsize) .- halo_pos
        # end

        # # transform velocities into halo frame of reference
        # if "VEL" in blocks && length(p.vel) > 0
            # p.vel .-= halo_vel
        # end

        # convert_units!(p, h, units)

        # # read particle mass if available in header mass array
        # if "MASS" ∉ keys(p)
            # p_mass = read_header_particle_mass(g.snapshot, ptype, units)
            # # ustrip removes any potential unit to be comparable with 0 (Int64)
            # ustrip(p_mass) > 0 && (p.mass = p_mass)
        # end

        # # add particles to galaxy
        # g[ptype] = p
    # end

    # return g
# end


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
