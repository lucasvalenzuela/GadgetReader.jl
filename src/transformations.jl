@doc raw"""
    rotate(p::Particles, rotmat::AbstractMatrix{<:Real})
    rotate!(p::Particles, rotmat::AbstractMatrix{<:Real})

Rotates the particles `p` by the ``3 \times 3`` rotation matrix `rotmat`.

The non-in-place version `rotate` creates a copy of the particles with only new pointers to the
rotated properties (generally "POS" and "VEL").
"""
function rotate(p::Particles, rotmat::AbstractMatrix{<:Real})
    @assert size(rotmat) == (3, 3)

    pc = copy(p)

    # only rotate existing quantities (currently only pos and vel)
    for prop in intersect(keys(pc), [:pos, :vel])
        pc[prop] = rotate(pc[prop], rotmat)
    end

    return pc
end

function LinearAlgebra.rotate!(p::Particles, rotmat::AbstractMatrix{<:Real})
    @assert size(rotmat) == (3, 3)

    # only rotate existing quantities (currently only pos and vel)
    for prop in intersect(keys(p), [:pos, :vel])
        rotate!(p[prop], rotmat)
    end

    return p
end

@doc raw"""
    rotate(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})
    rotate!(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})

Rotates the galaxy `g` by the ``3 \times 3`` rotation matrix `rotmat`.

The non-in-place version `rotate` creates a copy of the galaxy with only new pointers to the
rotated quantites (generally "POS" and "VEL") and copies of the particle Dict.
"""
function rotate(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})
    @assert size(rotmat) == (3, 3)

    gc = copy(g)

    # do for stars, dm, etc.
    for p in values(gc)
        # only rotate existing quantities (currently only pos and vel)
        for prop in intersect(keys(p), [:pos, :vel])
            p[prop] = rotate(p[prop], rotmat)
        end
    end

    return gc
end

function LinearAlgebra.rotate!(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})
    @assert size(rotmat) == (3, 3)
    # do for stars, dm, etc.
    for p in values(g)
        rotate!(p, rotmat)
    end

    return g
end


"""
    rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

Returns `vals` rotated by the rotation matrix `rotmat`. Works for any dimensions.
"""
rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real}) = rotmat * vals

function LinearAlgebra.rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    vals .= rotate(vals, rotmat)
    return vals
end


"""
    rotate_edgeon(p::Particles; axis_ratios::Bool=false, kwargs...)
    rotate_edgeon!(p::Particles; axis_ratios::Bool=false, kwargs...)

Rotates the particles edge-on using the given shape determination algorithm and `radius` (see [`rotation_matrix_edgeon`](@ref) for all keyword parameters).

Keyword parameters
- `axis_ratios::Bool`: `true` to return tuple `(p, q, s)` with particles and axis ratios
"""
function rotate_edgeon(p::Particles; axis_ratios::Bool=false, kwargs...)
    Q⁻¹, q, s = rotation_matrix_edgeon(p; kwargs...)

    pc = rotate(p, Q⁻¹)

    if axis_ratios
        return pc, q, s
    end

    return pc
end

function rotate_edgeon!(p::Particles; axis_ratios::Bool=false, kwargs...)
    Q⁻¹, q, s = rotation_matrix_edgeon(p; kwargs...)

    rotate!(p, Q⁻¹)

    if axis_ratios
        return p, q, s
    end

    return p
end


"""
    rotate_edgeon(g::AbstractGalaxy, ptype::Symbol=:stars; axis_ratios::Bool=false, kwargs...)
    rotate_edgeon!(g::AbstractGalaxy, ptype::Symbol=:stars; axis_ratios::Bool=false, kwargs...)

Rotates the galaxy's particles edge-on for the given particle type `ptype` (`:stars`, `:dm`, `:gas`, etc.)
using the given shape determination algorithm and `radius`
(see [`rotation_matrix_edgeon`](@ref) for all keyword parameters).

Keyword parameters
- `axis_ratios::Bool`: `true` to return tuple `(g, q, s)` with galaxy and axis ratios
"""
function rotate_edgeon(g::AbstractGalaxy, ptype::Symbol=:stars; axis_ratios::Bool=false, kwargs...)
    Q⁻¹, q, s = rotation_matrix_edgeon(g[ptype]; kwargs...)

    gc = rotate(g, Q⁻¹)

    if axis_ratios
        return gc, q, s
    end

    return gc
end

function rotate_edgeon!(g::AbstractGalaxy, ptype::Symbol=:stars; axis_ratios::Bool=false, kwargs...)
    Q⁻¹, q, s = rotation_matrix_edgeon(g[ptype]; kwargs...)

    rotate!(g, Q⁻¹)

    if axis_ratios
        return g, q, s
    end

    return g
end
