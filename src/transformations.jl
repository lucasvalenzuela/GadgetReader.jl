@doc raw"""
    rotate!(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})

Rotates the galaxy `g` by the ``3 \times 3`` rotation matrix `rotmat`.
"""
function rotate!(g::AbstractGalaxy, rotmat::AbstractMatrix{<:Real})
    @assert size(rotmat) == (3, 3)
    # do for stars, dm, etc.
    for ptype in keys(g)
        # only rotate existing quantities (currently only pos and vel)
        for prop in intersect(keys(g[ptype]), [:pos, :vel])
            rotate!(g[ptype][prop], rotmat)
        end
    end

    return g
end


"""
    rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

Returns `vals` rotated by the rotation matrix `rotmat`. Works for any dimensions.
"""
rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real}) = rotmat * vals

function rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    vals .= rotate(vals, rotmat)
    return vals
end
