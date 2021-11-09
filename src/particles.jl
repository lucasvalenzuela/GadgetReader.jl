##################################
### Particles
##################################

"""
    Particles(type::Symbol, properties::Dict{String})

Particles of a type (typically `:stars`, `:dm`, `:gas`, or `:bh`) with their properties in a `Dict`.
The `Dict` has `String` keys, which are set to all uppercase by default, except when there is at
least one uppercase character in the key. In that case the key is left as is.

The properties are accessible via `particles.prop` (case insensitive).

# Examples
```jldoctest
julia> p = GadgetReader.Particles(:stars, Dict("ID"=>[1, 2], "MASS"=>[1e11, 2e11]));

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
        # completely lowercase is converted to uppercase, otherwise is assumed to be correct
        sym_str = String(sym)
        if all(c -> islowercase(c), sym_str)
            key = uppercase(sym_str)
        else
            key = sym_str
        end
        return obj.properties[key]
    end
end

function Base.setproperty!(obj::Particles, sym::Symbol, val)
    if sym in fieldnames(Particles)
        setfield!(obj, sym, val)
    else
        # completely lowercase is converted to uppercase, otherwise is assumed to be correct
        sym_str = String(sym)
        if all(c -> islowercase(c), sym_str)
            key = uppercase(sym_str)
        else
            key = sym_str
        end
        obj.properties[key] = val
    end
end

Base.getindex(obj::Particles, str::String) = obj.properties[str]
Base.getindex(obj::Particles, sym::Symbol) = Base.getproperty(obj, sym)
Base.setindex!(obj::Particles, val, str::String) = Base.setindex!(obj.properties, val, str)
Base.setindex!(obj::Particles, val, sym::Symbol) = Base.setproperty!(obj, sym, val)
function Base.keys(obj::Particles)
    [all(c -> isuppercase(c), key) ? Symbol(lowercase(key)) : Symbol(key) for key in keys(obj.properties)]
end
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
    props = join(keys(obj.properties) .|> String |> sort, " ")
    print(io, props)
end

"""
    particle_type_id(config::GadgetConfig, type::Symbol)

Returns the Gadget particle type from a particle `Symbol`, according to the GADGET config.
"""
particle_type_id(config::GadgetConfig, type::Symbol) = config.particle_types[type]
