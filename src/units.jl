function convert_units!(p::Particles, h::SnapshotHeader, units::Symbol=:full)
    if units === :physical
        for (key, vals) in pairs(p)
            convert_units_physical!(vals, key, h)
        end
    elseif units === :full
        for (key, vals) in pairs(p)
            p[key] = convert_units_full(vals, key, h)
        end
    end

    props = keys(p)
    if :zs in props && (:im in props || :mass in props)
        p.zs = convert_units_solar_metallicity(p.zs, :im in props ? p.iM : p.mass)
    end

    return p
end

for (type, excl) in [("physical", ""), ("physical", "!"), ("full", "")]
    quote
        function $(Symbol("convert_units_", type, excl))(
            vals::AbstractArray{<:Real},
            prop::Symbol,
            h::SnapshotHeader,
        )
            if prop === :pos
                $(Symbol("convert_units_", type, "_pos", excl))(vals, h)
            elseif prop === :vel
                $(Symbol("convert_units_", type, "_vel", excl))(vals, h)
            elseif prop === :age
                $(Symbol("convert_units_", type, "_age", excl))(vals, h)
            elseif prop === :temp
                $(Symbol("convert_units_", type, "_temp", excl))(vals, h)
            elseif prop in [:mass, :zs, :im]
                $(Symbol("convert_units_", type, "_mass", excl))(vals, h)
            else
                vals
            end
        end
    end |> eval
end

for (key, factor, unit) in [
    ("pos", :(1 / (h.h0 * (h.z + 1))), :(u"kpc")),
    ("vel", :(1 / sqrt(h.z + 1)), :(u"km/s")),
    ("temp", :(1), :(u"K")),
    ("mass", :(1e10 / h.h0), :(u"Msun")),
]
    quote
        function $(Symbol("convert_units_physical_", key))(val::T, h::SnapshotHeader) where {T<:Real}
            val * convert(T, $factor)
        end
        function $(Symbol("convert_units_full_", key))(val::Real, h::SnapshotHeader)
            $(Symbol("convert_units_physical_", key))(val, h) * $unit
        end
        function $(Symbol("convert_units_physical_", key))(
            vals::AbstractArray{T},
            h::SnapshotHeader,
        ) where {T<:Real}
            vals .* convert(T, $factor)
        end
        function $(Symbol("convert_units_physical_", key, "!"))(
            vals::AbstractArray{T},
            h::SnapshotHeader,
        ) where {T<:Real}
            vals .*= convert(T, $factor)
        end
        function $(Symbol("convert_units_full_", key))(
            vals::AbstractArray{T},
            h::SnapshotHeader,
        ) where {T<:Real}
            vals .* (convert(T, $factor) * $unit)
        end
    end |> eval
end

function convert_units_physical_age(val::Real, h::SnapshotHeader)
    convert_units_full_age(val, h) |> ustrip
end
function convert_units_full_age(val::T, h::SnapshotHeader) where {T<:Real}
    c = cosmology(h=h.h0, OmegaM=h.omega_0)
    lookback_time(c, 1 / val - 1) - lookback_time(c, h.z) # in Gyr
end
function convert_units_physical_age(vals::AbstractArray{<:Real}, h::SnapshotHeader)
    convert_units_full_age(vals, h) .|> ustrip
end
function convert_units_physical_age!(vals::AbstractArray{<:Real}, h::SnapshotHeader)
    vals .= convert_units_physical_age(vals, h)
end
function convert_units_full_age(vals::AbstractArray{T}, h::SnapshotHeader) where {T<:Real}
    c = cosmology(h=h.h0, OmegaM=h.omega_0)
    lookback_time.((c,), 1 ./ vals .- 1) .- lookback_time(c, h.z) # in Gyr
end


function convert_units_solar_metallicity(zs::AbstractMatrix{<:Number}, mass::AbstractVector{<:Number})
    dropdims(sum(zs[2:11, :]; dims=1); dims=1) ./ (mass .- dropdims(sum(zs; dims=1); dims=1)) ./ 0.0134 # in solar metallicities
end
