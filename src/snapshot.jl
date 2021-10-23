"""
    Snapshot(snapbase, subbase)

Returns a Snapshot from a `snapbase` and a `subbase`, which can also be `nothing`.
Both can be files or the base part of files ending with `.0`, `.1`, etc.
"""
struct Snapshot{T}
    snapbase::Union{T,Nothing}
    subbase::Union{T,Nothing}
end

"""
    Snapshot(; snapbase=nothing, subbase=nothing)

Convenience method for creating a `Snapshot` where one of the base files does not exist.
"""
Snapshot(; snapbase=nothing, subbase=nothing) = Snapshot(snapbase, subbase)

"""
    Snapshot(path_box, snap::Integer; snapbase=true, subbase=true)

Convenience method for creating a `Snapshot` from the path to a simulated box and a snapshot number.
Setting `snapbase` or `subbase` to `false` sets the respective file base to `nothing`.

The formats used are the following, for `snapbase` and `subbase`, respectively:
- `snapdir_XXX/snap_XXX`
- `groups_XXX/sub_XXX`
"""
function Snapshot(path_box, snap::Integer; snapbase::Bool=true, subbase::Bool=true)
    _snapbase = snapbase ? get_snapbase(path_box, snap) : nothing
    _subbase = subbase ? get_subbase(path_box, snap) : nothing
    Snapshot(_snapbase, _subbase)
end

get_snapbase(path_box, snap::String) = joinpath(path_box, "snapdir_$snap", "snap_$snap")
get_snapbase(path_box, snap::Integer) = get_snapbase(path_box, lpad(snap, 3, '0'))
get_subbase(path_box, snap::String) = joinpath(path_box, "groups_$snap", "sub_$snap")
get_subbase(path_box, snap::Integer) = get_subbase(path_box, lpad(snap, 3, '0'))

function Base.show(io::IO, ::MIME"text/plain", obj::Snapshot)
    printstyled(io, Snapshot; bold=true)
    println()
    isnothing(obj.snapbase) || println(io, "Snap: $(obj.snapbase)")
    isnothing(obj.subbase) || print(io, "Sub:  $(obj.subbase)")
end
