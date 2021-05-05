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
Base.setindex!(obj::Particles, val, sym::Symbol) = Base.setproperty!(obj, sym, val)
Base.keys(obj::Particles) = keys(obj.properties) .|> lowercase .|> Symbol
Base.values(obj::Particles) = values(obj.properties)
Base.propertynames(obj::Particles) = [:type; keys(obj)]

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
Will try to read "MSUB" to extract the  if `get_id` is `true`.

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
Base.keys(obj::Galaxy) = keys(obj.particles) |> collect
Base.values(obj::Galaxy) = values(obj.particles)
Base.propertynames(obj::Galaxy) = [[:snapshot, :isub, :subid]; keys(obj)]

function Base.show(io::IO, ::MIME"text/plain", obj::Galaxy)
    printstyled(io, "Galaxy"; bold=true)
    try
        z = round(read_redshift(obj.snapshot); digits=2)
        println(io, " at z=$z")
    catch
        println()
    end
    println(io, " isub $(obj.isub)")
    for key in sort(keys(obj))
        show(io, "text/plain", obj[key])
    end
end
