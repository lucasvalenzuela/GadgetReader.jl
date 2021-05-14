##################################
### Particles
##################################

"""
    Particles(type::Symbol, properties::Dict{String})

Particles of a type (typically `:stars`, `:dm`, `:gas`, or `:bh`) with their properties in a `Dict`.
The `Dict` has `String` keys, which are set to all uppercase by default, with exception
of `iM` and `Zs`.

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
        return obj.properties[sym |> String |> uppercase |> _particle_property_exceptions]
    end
end

function Base.setproperty!(obj::Particles, sym::Symbol, val)
    if sym in fieldnames(Particles)
        setfield!(obj, sym, val)
    else
        obj.properties[sym |> String |> uppercase |> _particle_property_exceptions] = val
    end
end

Base.getindex(obj::Particles, str::String) = obj.properties[str]
Base.getindex(obj::Particles, sym::Symbol) = Base.getproperty(obj, sym)
Base.setindex!(obj::Particles, val, str::String) = Base.setindex!(obj.properties, val, str)
Base.setindex!(obj::Particles, val, sym::Symbol) = Base.setproperty!(obj, sym, val)
Base.keys(obj::Particles) = keys(obj.properties) .|> lowercase .|> Symbol
Base.haskey(obj::Particles, key) = haskey(obj.properties, key)
Base.haskey(obj::Particles, key::Symbol) = key in keys(obj)
Base.values(obj::Particles) = values(obj.properties)
Base.propertynames(obj::Particles) = [:type; keys(obj)]
Base.copy(obj::Particles) = Particles(obj.type, copy(obj.properties))

function Base.show(io::IO, ::MIME"text/plain", obj::Particles)
    printstyled(io, String(obj.type); bold=true)
    len = 0
    for val in values(obj)
        if typeof(val) <: AbstractArray && size(val, ndims(val)) > len
            len = size(val, ndims(val))
        end
    end
    println(io, ": $len Particles")
    print(io, " ")
    props = join(keys(obj.properties) .|> _particle_property_exceptions |> sort, " ")
    println(io, props)
end

function _particle_property_exceptions(key::String)
    if key == "IM"
        return "iM"
    elseif key == "ZS"
        return "Zs"
    else
        return key
    end
end


const _particle_type_id = Dict(:gas => 0, :dm => 1, :stars => 4, :bh => 5)

"""
    particle_type_id(type::Symbol)

Returns the Gadget particle type from a particle `Symbol`
(currently `:gas`, `:dm`, `:stars`, and `:bh`).
"""
particle_type_id(type::Symbol) = _particle_type_id[type]



##################################
### Galaxy
##################################

abstract type AbstractGalaxy end

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
struct Galaxy <: AbstractGalaxy
    snapshot::Snapshot
    isub::Integer
    subid::Union{HaloID,Nothing}
    particles::Dict{Symbol,Particles}
end

"""
    Galaxy(snapshot::Snapshot, isub[, get_id=true])

Convenience method for creating a `Galaxy` with only a `snapshot` and a subfind id.
Will try to read "MSUB" to extract the  if `get_id` is `true`.

[`Particles`](@ref) can still be added after initializing the galaxy by calling
for example `galaxy.stars = [...]`.
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

"""
    GalaxyGroup(snapshot::Snapshot,
           igroup::Integer,
           groupid::Union{HaloID,Nothing},
           particles::Dict{Symbol,Particles})

Group of a given `snapshot` ([`Snapshot`](@ref)) with the zero-based group index `igroup`
with [`Particles`](@ref).
Including `groupid` (containing information on the exact position in the subfind files) makes
accessing properties of the galaxy faster than when only having `igroup`.

The particles are accessible via `group.particletype` or `group[:particletype]`
(see [`Particles`](@ref) for typical types).

# Examples
```jldoctest
julia> p = GadgetGalaxies.Particles(:stars, Dict("ID"=>[1, 2], "MASS"=>[1e11, 2e11]));

julia> g = GalaxyGroup(Snapshot("box", 13), 164, nothing, Dict(:stars=>p));

julia> g.stars === p
true
```
"""
struct GalaxyGroup <: AbstractGalaxy
    snapshot::Snapshot
    igroup::Integer
    groupid::Union{HaloID,Nothing}
    particles::Dict{Symbol,Particles}
end

"""
    GalaxyGroup(snapshot::Snapshot, igroup[, get_id=true])

Convenience method for creating a `GalaxyGroup` with only a `snapshot` and a subfind id.
Will try to read "MTOT" to extract the  if `get_id` is `true`.

[`Particles`](@ref) can still be added after initializing the group by calling for example
`group.stars = [...]`.
"""
function GalaxyGroup(snapshot::Snapshot, igroup::Integer, get_id::Bool)
    if get_id
        subbase = string(snapshot.subbase)
        subfile = GadgetIO.select_file(subbase, 0)

        n_files = read_header(subfile).num_files
        _, groupid = read_halo_prop_and_id(subbase, igroup, "MTOT", n_files; verbose=false)
    else
        groupid = nothing
    end

    return GalaxyGroup(snapshot, igroup, groupid, Dict{Symbol,Particles}())
end

GalaxyGroup(snapshot::Snapshot, igroup::Integer) = GalaxyGroup(snapshot, igroup, true)


function Base.getproperty(obj::T, sym::Symbol) where {T<:AbstractGalaxy}
    if sym in fieldnames(T)
        return getfield(obj, sym)
    else
        return obj.particles[sym]
    end
end

function Base.setproperty!(obj::T, sym::Symbol, val) where {T<:AbstractGalaxy}
    if sym in fieldnames(T)
        setfield!(obj, sym, val)
    else
        obj.particles[sym] = val
    end
end

Base.getindex(obj::AbstractGalaxy, sym::Symbol) = Base.getproperty(obj, sym)
Base.setindex!(obj::AbstractGalaxy, val, sym::Symbol) = Base.setproperty!(obj, sym, val)
Base.keys(obj::AbstractGalaxy) = keys(obj.particles) |> collect
Base.haskey(obj::AbstractGalaxy, key) = haskey(obj.particles, key)
Base.values(obj::AbstractGalaxy) = values(obj.particles)
Base.propertynames(obj::Galaxy) = [[:snapshot, :isub, :subid]; keys(obj)]
Base.propertynames(obj::GalaxyGroup) = [[:snapshot, :igroup, :groupid]; keys(obj)]

function Base.copy(obj::T) where {T<:AbstractGalaxy}
    particles = Dict(key => copy(p) for (key, p) in pairs(obj.particles))
    T(obj.snapshot, getind(obj), getid(obj), particles)
end

function Base.show(io::IO, ::MIME"text/plain", obj::AbstractGalaxy)
    printstyled(io, typeof(obj); bold=true)
    try
        z = round(read_redshift(obj.snapshot); digits=2)
        println(io, " at z=$z")
    catch
        println()
    end
    println(io, " $(getindname(obj)) $(getind(obj))")
    for key in sort(keys(obj))
        show(io, "text/plain", obj[key])
    end
end


getindname(g::Galaxy) = "isub"
getindname(g::GalaxyGroup) = "igroup"

getind(g::Galaxy) = g.isub
getind(g::GalaxyGroup) = g.igroup

getid(g::Galaxy) = g.subid
getid(g::GalaxyGroup) = g.groupid

get_subfind_type(g::Galaxy) = 2
get_subfind_type(g::GalaxyGroup) = 1
