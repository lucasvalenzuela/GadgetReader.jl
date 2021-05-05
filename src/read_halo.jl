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
- `radius::Union{Real,Nothing}=nothing`: if `nothing`, the half-mass radius "RHMS"
- `rad_scale::Real=1`: read particles within `rad_scale * radius` of halo position
- `units::Symbol=:full`: `:full` for units as `Unitful` quantities, `:physical` for values converted
  to physical units (kpc, km/s, solar metallicities etc.), `:sim` for values in simulation units
- `use_keys::Bool=true`: if Peano-Hilbert keys should be used if possible
- `verbose::Bool=false`: verbose output
- `props::Tuple=(
       (:gas, ["POS", "VEL", "MASS", "TEMP", "SFR", "Zs"]),
       (:dm, ["POS", "VEL"]),
       (:stars, ["POS", "VEL", "MASS", "AGE", "Zs", "iM"]),
   )`: particle types (see [`Particles`](@ref)) and properties to be read
"""
function read_halo!(
    g::Galaxy;
    radius::Union{Real,Nothing}=nothing,
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

    # Provides properties `h0`, `z`, `time`, `omega_0`, `omega_l`, `num_files`
    h = read_header(GadgetIO.select_file(subbase, 0))

    # get global halo properties
    halo_pos = read_halo_prop(subbase, g.subid, "SPOS"; verbose)
    halo_vel = read_halo_prop(subbase, g.subid, "SVEL"; verbose)

    # handle particles for different types
    Threads.@threads for (ptype, blocks) in props
        ptype_id = particle_type_id(ptype)

        prop_data = read_particles_in_halo(
            snapbase,
            blocks,
            subbase,
            g.subid;
            parttype=ptype_id,
            halo_type=2,
            rad_scale,
            use_keys,
            verbose,
        )
        p = Particles(ptype, prop_data)

        # shift across box borders and transform into halo frame of reference
        if "POS" in blocks
            p.pos .= shift_across_box_border.(p.pos, halo_pos, h.boxsize, 1 // 2 * h.boxsize) .- halo_pos
        end

        # transform velocities into halo frame of reference
        if "VEL" in blocks
            p.vel .-= halo_vel
        end

        convert_units!(p, h, units)

        # add particles to galaxy
        g[ptype] = p
    end

    return g
end


"""
    read_redshift(snapshot::Snapshot)

Returns the redshift `z` from the first snapfile's header.
"""
function read_redshift(snapshot::Snapshot)
    h = read_subfind_header(GadgetIO.select_file(snapshot.subbase |> string, 0))
    return h.z
end
