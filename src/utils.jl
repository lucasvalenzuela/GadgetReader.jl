skipnan(x) = Iterators.filter(!isnan, x)

function get_dims(perspective::Symbol)
    dims = if perspective === :edgeon
        [1, 3]
    elseif perspective === :sideon
        [2, 3]
    elseif perspective === :faceon
        [1, 2]
    else
        error("`perspective` is $perspective, needs to be one of `:edgeon`, `:sideon`, `:faceon`")
    end
end

function LinearAlgebra.cross(a::AbstractMatrix, b::AbstractMatrix)
    @assert size(a) == size(b)
    @assert size(a, 1) == 3

    c = Array{typeof(a[1] * b[1])}(undef, 3, size(a, 2))
    @views @. c[1, :] = a[2, :] * b[3, :] - a[3, :] * b[2, :]
    @views @. c[2, :] = a[3, :] * b[1, :] - a[1, :] * b[3, :]
    @views @. c[3, :] = a[1, :] * b[2, :] - a[2, :] * b[1, :]

    return c
end
