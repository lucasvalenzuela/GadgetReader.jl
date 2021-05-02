##################################
### Particles
##################################

"""
    Particles(type::Symbol, properties::Dict{String})

Particles of a type (typically `:stars`, `:dm`, `:gas`, or `:bh`) with their properties in a [`Dict`](@ref).
The [`Dict`](@ref) has [`String`](@ref) keys, which are set to all uppercase by default.

The properties are accessible via `particles.prop` (case insensitive).

# Examples
```jldoctest
julia> p = GadgetGalaxies.Particles(:stars, Dict("ID"=>[1, 2], "MASS"=>[1e11, 2e11]));

julia> p.id
2-element Vector{Int64}:
 1
 2

julia> p.pos = [1 2; 3 4; 5 6]
3×2 Matrix{Int64}:
 1  2
 3  4
 5  6

julia> p.properties["POS"] === p.pos
true
```
"""
struct Particles
    type::Symbol
    properties::Dict{String,Any}
end

function Base.getproperty(obj::Particles, sym::Symbol)
    if sym in fieldnames(Particles)
        return getfield(obj, sym)
    else
        return obj.properties[sym |> String |> uppercase]
    end
end

function Base.setproperty!(obj::Particles, sym::Symbol, val)
    if sym in fieldnames(Particles)
        setfield!(obj, sym, val)
    else
        obj.properties[sym |> String |> uppercase] = val
    end
end

Base.getindex(obj::Particles, sym::Symbol) = Base.getproperty(obj, sym)
Base.setindex!(obj::Particles, val, sym::Symbol) = Base.setproperty!(obj, sym, val)

"""
    particleproperties(particles::Particles)

Returns [`Vector`](@ref) of particle properties saved in `particles` (e.g. `[:id, :pos, :vel]`).
"""
particleproperties(obj::Particles) = keys(obj.properties) .|> lowercase .|> Symbol

function Base.propertynames(obj::Particles)
    return [:type; particleproperties(obj)]
end



##################################
### Galaxy
##################################

"""
    Galaxy(snapshot::Snapshot,
           isub::Integer,
           subid::Union{HaloID,Nothing},
           particles::Dict{Symbol,Particles})

Galaxy of a given `snapshot` ([`Snapshot`](@ref)) with the zero-based subhalo index `isub`
with [`Particles`](@ref).
Including `subid` (containing information on the exact position in the subfind files) makes
accessing properties of the galaxy faster than when only having `isub`.

The particles are accessible via `galaxy.particletype` or `galaxy[:particletype]`
(see [`Particles`](@ref) for typical types).

# Examples
```jldoctest
julia> p = GadgetGalaxies.Particles(:stars, Dict("ID"=>[1, 2], "MASS"=>[1e11, 2e11]));

julia> g = Galaxy(Snapshot("box", 13), 1532, nothing, Dict(:stars=>p));

julia> g.stars === p
true
```
"""
struct Galaxy
    snapshot::Snapshot
    isub::Integer
    subid::Union{HaloID,Nothing}
    particles::Dict{Symbol,Particles}
end

"""
    Galaxy(snapshot::Snapshot, isub[, get_id=true])

Convenience method for creating a `Galaxy` with only a `snapshot` and a subfind id.
Will try to read "MSUB" to extract the [`HaloID`](@ref) if `get_id` is `true`.

[`Particles`](@ref) can still be added after initializing the halo by calling for example `galaxy.stars = [...]`.
"""
function Galaxy(snapshot::Snapshot, isub::Integer, get_id::Bool)
    if get_id
        subbase = string(snapshot.subbase)
        subfile = GadgetIO.select_file(subbase, 0)

        n_files = read_header(subfile).num_files
        _, subid = read_halo_prop_and_id(subbase, isub, "MSUB", n_files; verbose=false)
    else
        subid = nothing
    end

    return Galaxy(snapshot, isub, subid, Dict{Symbol,Particles}())
end

Galaxy(snapshot::Snapshot, isub::Integer) = Galaxy(snapshot, isub, true)

function Base.getproperty(obj::Galaxy, sym::Symbol)
    if sym in fieldnames(Galaxy)
        return getfield(obj, sym)
    else
        return obj.particles[sym]
    end
end

function Base.setproperty!(obj::Galaxy, sym::Symbol, val)
    if sym in fieldnames(Galaxy)
        setfield!(obj, sym, val)
    else
        obj.particles[sym] = val
    end
end

Base.getindex(obj::Galaxy, sym::Symbol) = Base.getproperty(obj, sym)
Base.setindex!(obj::Galaxy, val, sym::Symbol) = Base.setproperty!(obj, sym, val)

"""
    particletypes(galaxy::Galaxy)

Returns a `Vector` of particle types saved in `galaxy` (e.g. `[:stars, :dm]`).
"""
particletypes(obj::Galaxy) = keys(obj.particles) |> collect

Base.propertynames(obj::Galaxy) = [[:snapshot, :isub, :subid]; particletypes(obj)]
