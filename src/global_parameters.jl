"""
    half_mass_radius(
        pos::AbstractMatrix{<:Number},
        mass::Union{AbstractVector{<:Number},Nothing},
        [r²::AbstractVector{<:Number}];
        rmax::Union{Number,Nothing}=nothing,
    )

Returns the half-mass radius of the particles within radius `rmax`, or of all particles if `rmax`
is `nothing`.
"""
function half_mass_radius(
    pos::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number},
    r²::AbstractVector{<:Number};
    rmax::Union{Number,Nothing}=nothing,
)
    # check dimensions
    @assert size(pos, 1) == 3
    @assert length(mass) == size(pos, 2) == length(r²)

    if isnothing(rmax)
        mask = (:) # all particles
    else
        mask = r² .≤ rmax^2
    end

    mass_half = @views 1 // 2 * sum(mass[mask])

    # mass sorted by distance to center
    ind = @views sortperm(r²[mask])
    mass_sorted = mass[mask][ind]

    # traverse masses outwardly until half of the total mass is reached
    mass_cumul = zero(eltype(mass))
    i = 1
    @inbounds while mass_cumul < mass_half
        mass_cumul += mass_sorted[i]
        i += 1
    end

    r²_masked = @views r²[mask]
    return 1 // 2 * (sqrt(r²_masked[ind[i - 1]]) + sqrt(r²_masked[ind[i - 2]]))
end

function half_mass_radius(
    pos::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number};
    rmax::Union{Number,Nothing}=nothing,
)
    r² = r²_sphere(pos)
    half_mass_radius(pos, mass, r²; rmax)
end

"""
    half_mass_radius(pos::AbstractMatrix{<:Number}; rmax::Union{Number,Nothing}=nothing)

Returns the half-mass radius of the particles within radius `rmax` assuming equal masses of all
particles, or of all particles if `rmax` is `nothing`.
"""
function half_mass_radius(pos::AbstractMatrix{<:Number}; rmax::Union{Number,Nothing}=nothing)
    # check dimensions
    @assert size(pos, 1) == 3

    r² = r²_sphere(pos)

    if isnothing(rmax)
        return quantile(sqrt.(r²), 0.5)
    else
        mask = r² .≤ rmax^2
        return @views quantile(sqrt.(r²[mask]), 0.5)
    end
end

half_mass_radius(pos::AbstractMatrix{<:Number}, mass::Nothing; kwargs...) = half_mass_radius(pos, kwargs...)

"""
    half_mass_radius(g::AbstractGalaxy, ptype::Symbol; kwargs...)
    half_mass_radius(p::Particles; kwargs...)

Returns the particles' half-mass radius, depending on if a mass vector is available or not.
If the mass is not available, equal particle masses are assumed.
"""
function half_mass_radius(p::Particles; kwargs...)
    if haskey(p, :mass) && p.mass isa AbstractVector
        half_mass_radius(p.pos, p.mass; kwargs...)
    else
        half_mass_radius(p.pos; kwargs...)
    end
end

half_mass_radius(g::AbstractGalaxy, ptype::Symbol; kwargs...) = half_mass_radius(g[ptype]; kwargs...)


"""
    half_mass_radius_2D(
        pos::AbstractMatrix{<:Number},
        [mass::Union{AbstractVector{<:Number},Nothing}];
        rmax::Union{Number,Nothing}=nothing,
        perspective::Symbol=:edgeon,
    )

Returns the 2D half-mass radius of the particles within radius `rmax` when viewed from a given
`perspective` (`:edgeon`, `:sideon`, `:faceon`), or of all particles if `rmax` is `nothing`.
"""
function half_mass_radius_2D(
    pos::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number};
    rmax::Union{Number,Nothing}=nothing,
    perspective::Symbol=:edgeon,
)
    # check dimensions
    @assert size(pos, 1) == 3

    dims = get_dims(perspective)
    r² = @views r²_circle(pos[dims, :])

    return half_mass_radius(pos, mass, r²; rmax)
end

function half_mass_radius_2D(
    pos::AbstractMatrix{<:Number};
    rmax::Union{Number,Nothing}=nothing,
    perspective::Symbol=:edgeon,
)
    # check dimensions
    @assert size(pos, 1) == 3

    dims = get_dims(perspective)
    r² = @views r²_circle(pos[dims, :])

    if isnothing(rmax)
        return quantile(sqrt.(r²), 0.5)
    else
        mask = r² .≤ rmax^2
        return @views quantile(sqrt.(r²[mask]), 0.5)
    end
end

function half_mass_radius_2D(pos::AbstractMatrix{<:Number}, mass::Nothing; kwargs...)
    half_mass_radius_2D(pos; kwargs...)
end

"""
    half_mass_radius_2D(g::AbstractGalaxy, ptype::Symbol; kwargs...)
    half_mass_radius_2D(p::Particles; kwargs...)

Returns the particles' 2D half-mass radius, depending on if a mass vector is available or not.
If the mass is not available, equal particle masses are assumed.
"""
function half_mass_radius_2D(p::Particles; kwargs...)
    if haskey(p, :mass) && p.mass isa AbstractVector
        half_mass_radius_2D(p.pos, p.mass; kwargs...)
    else
        half_mass_radius_2D(p.pos; kwargs...)
    end
end

half_mass_radius_2D(g::AbstractGalaxy, ptype::Symbol; kwargs...) = half_mass_radius_2D(g[ptype]; kwargs...)
