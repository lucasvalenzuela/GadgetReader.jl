"""
    convert_units!(p::Particles, h, units::Symbol=:full)

Converts the units of the particles properties (see [`Particles`](@ref)) from a snapshot or
subfind header in-place.
Use `:full` for units as `Unitful` quantities, `:physical` for values converted
to physical units (kpc, km/s, solar metallicities etc.), or `:sim` for values in simulation units.
The individual metallicities in "Zs" are converted to a one-dimensional `Vector` of values in solar
metallicities in any case.
"""
function convert_units!(p::Particles, h, units::Symbol=:full)
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
    if :zs in props && length(p.zs) > 0 && (:im in props || :mass in props)
        p.zs = convert_units_solar_metallicity(p.zs, :im in props ? p.iM : p.mass)
    end

    return p
end

function convert_units_subfind_prop(
    val::Union{Real,AbstractArray{<:Real}},
    prop::AbstractString,
    h,
    units::Symbol=:full;
    verbose::Bool=false,
)
    # TODO: SPIN (units: angular momentum), DSUB (vel. dispersion units?),
    # SLUM, SLAT, SLOB, DUST, SZ (units?), SSFR (units really M⊙/yr?)
    if prop[1] === 'R' || prop in ["GPOS", "BGPO", "BGRA", "SPOS", "SCM", "SHMR"]
        return convert_units_pos(val, h, units)
    elseif prop[1] === 'V' || prop in ["SVEL"]
        return convert_units_vel(val, h, units)
    elseif (prop[1] === 'M' && prop != "MBID") || prop in ["BGMA", "SMST"]
        return convert_units_mass(val, h, units)
    elseif prop in ["SAGE"]
        return convert_units_age(val, h, units)
        # elseif prop in ["TGAS"] # TODO: unclear what KeV unit is
        # return convert_units_temp(val, h, units)
        # elseif prop in ["LGAS"] # TODO: energy unit conversion (10^44 erg/s)
        # return convert_units_temp(val, h, units)
    elseif prop in ["SSFR"] && units === :full
        return val * u"Msun" / Unitful.yr
    else
        verbose && @warn "The quantity $prop is returned in simulation units"
        return val
    end
end


# define generated functions for static linting
function convert_units_full end
function convert_units_physical end
function convert_units_physical! end

function simulation_units_pos end
function simulation_units_pos! end
function convert_units_pos end
function convert_units_full_pos end
function convert_units_physical_pos end
function convert_units_physical_pos! end
function simulation_units_vel end
function simulation_units_vel! end
function convert_units_vel end
function convert_units_full_vel end
function convert_units_physical_vel end
function convert_units_physical_vel! end
function simulation_units_temp end
function simulation_units_temp! end
function convert_units_temp end
function convert_units_full_temp end
function convert_units_physical_temp end
function convert_units_physical_temp! end
function simulation_units_mass end
function simulation_units_mass! end
function convert_units_mass end
function convert_units_full_mass end
function convert_units_physical_mass end
function convert_units_physical_mass! end

for (type, excl) in [("physical", ""), ("physical", "!"), ("full", "")]
    quote
        function $(Symbol("convert_units_", type, excl))(vals::AbstractArray{<:Real}, prop::Symbol, h)
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
@doc """
    convert_units_physical(vals::AbstractArray{<:Real}, prop::Symbol, h)
    convert_units_physical!(vals::AbstractArray{<:Real}, prop::Symbol, h)
    convert_units_full(vals::AbstractArray{<:Real}, prop::Symbol, h)

Converts simulation values to the respective physical values, depending on `prop`,
which can take any of the following values: `:pos`, `:vel`, `:temp`, `:mass`, `:age`.

Full returns values in `Unitful` quantities, whereas physical returns the physical value without unit.
"""
convert_units_physical,
#convert_units_physical!,
convert_units_full

for (key, factor, unit, plural, eq) in [
    ("pos", :(1 / (h.h0 * (h.z + 1))), :(u"kpc"), "positions", raw"x \to x/(h_0 (z+1)))"),
    ("vel", :(1 / sqrt(h.z + 1)), :(u"km/s"), "velocities", raw"v \to v/\sqrt{z+1})"),
    ("temp", :(1), :(u"K"), "temperatures", raw"T \to T"),
    ("mass", :(1e10 / h.h0), :(u"Msun"), "masses", raw"m \to m \times 10^{10} / h_0"),
]
    quote
        function $(Symbol("convert_units_physical_", key))(val::T, h) where {T<:Real}
            val * convert(T, $factor)
        end
        function $(Symbol("convert_units_full_", key))(val::Real, h)
            $(Symbol("convert_units_physical_", key))(val, h) * $unit
        end
        function $(Symbol("convert_units_physical_", key))(vals::AbstractArray{T}, h) where {T<:Real}
            vals .* convert(T, $factor)
        end
        function $(Symbol("convert_units_physical_", key, "!"))(vals::AbstractArray{T}, h) where {T<:Real}
            vals .*= convert(T, $factor)
        end
        function $(Symbol("convert_units_full_", key))(vals::AbstractArray{T}, h) where {T<:Real}
            vals .* (convert(T, $factor) * $unit)
        end
        function $(Symbol("convert_units_", key))(
            val::Union{Real,AbstractArray{<:Real}},
            h,
            units::Symbol=:full,
        )
            if units === :physical
                $(Symbol("convert_units_physical_", key))(val, h)
            elseif units === :full
                $(Symbol("convert_units_full_", key))(val, h)
            else
                val
            end
        end
        function $(Symbol("simulation_units_", key))(val::Quantity{T}, h) where {T}
            ustrip($unit, val) / convert(T, $factor)
        end
        function $(Symbol("simulation_units_", key))(val::T, h) where {T<:Real}
            val / convert(T, $factor)
        end
        function $(Symbol("simulation_units_", key))(vals::AbstractArray{<:Quantity{T}}, h) where {T}
            ustrip.($unit, vals) ./ convert(T, $factor)
        end
        function $(Symbol("simulation_units_", key))(vals::AbstractArray{T}, h) where {T<:Real}
            vals ./ convert(T, $factor)
        end
        function $(Symbol("simulation_units_", key, "!"))(vals::AbstractArray{T}, h) where {T<:Real}
            vals ./= convert(T, $factor)
        end
    end |> eval

    quote
        @doc """
             convert_units_$($key)(val::Union{Real,AbstractArray{<:Real}}, h, units::Symbol=:full)
             convert_units_physical_$($key)(val::Real, h)
             convert_units_physical_$($key)(vals::AbstractArray{<:Real}, h)
             convert_units_physical_$($key)!(vals::AbstractArray{<:Real}, h)
             convert_units_full_$($key)(val::Real, h)
             convert_units_full_$($key)(vals::AbstractArray{<:Real}, h)

         Converts $($plural) via ``$($eq)`` according to the cosmology defined by the snapshot header.
         `:full` returns values in `Unitful` quantities, whereas `:physical` returns the physical value
         in $(eval($unit)). `:sim` simply returns `val`.
         """
        $(Symbol("convert_units_", key))#,
        # $(Symbol("convert_units_physical_", key)),
        # $(Symbol("convert_units_physical_", key, "!")),
        # $(Symbol("convert_units_full_", key))

        @doc """
             simulation_units_$($key)(val::Number, h)
             simulation_units_$($key)(vals::AbstractArray{<:Number}, h)
             simulation_units_$($key)!(vals::AbstractArray{<:Real}, h)

         Converts $($plural) via ``$(replace($eq, "\\times " => "/ (") * (occursin("times", $eq) ? ")" : ""))``
         according to the cosmology defined by the snapshot header from physical to simulation units.
         Unitless parameters need to be given in $(eval($unit)).
         """
        $(Symbol("simulation_units_", key))#,
        # $(Symbol("simulation_units_", key, "!"))
    end |> eval
end

@doc raw"""
    convert_units_age(val::Union{Real,AbstractArray{<:Real}}, h, units::Symbol=:full)
    convert_units_physical_age(val::Real, h)
    convert_units_physical_age(vals::AbstractArray{<:Real}, h)
    convert_units_physical_age!(vals::AbstractArray{<:Real}, h)
    convert_units_full_age(val::Real, h)
    convert_units_full_age(vals::AbstractArray{<:Real}, h)

Converts ages from the scale factor ``a`` via the lookback time according to the cosmology
defined by the snapshot header.
Full returns values in `Unitful` quantities, whereas physical returns the physical value in Gyr.
"""
function convert_units_age(val::Union{Real,AbstractArray{<:Real}}, h, units::Symbol=:full)
    if units === :physical
        convert_units_physical_age(val, h)
    elseif units === :full
        convert_units_full_age(val, h)
    else
        val
    end
end
function convert_units_physical_age(val::Real, h)
    convert_units_full_age(val, h) |> ustrip
end
function convert_units_full_age(val::T, h) where {T<:Real}
    c = cosmology(h=h.h0, OmegaM=h.omega_0)
    lookback_time(c, 1 / val - 1) - lookback_time(c, h.z) # in Gyr
end
function convert_units_physical_age(vals::AbstractArray{<:Real}, h)
    convert_units_full_age(vals, h) .|> ustrip
end
function convert_units_physical_age!(vals::AbstractArray{<:Real}, h)
    vals .= convert_units_physical_age(vals, h)
end
function convert_units_full_age(vals::AbstractArray{T}, h) where {T<:Real}
    c = cosmology(h=h.h0, OmegaM=h.omega_0)
    lookback_time.((c,), 1 ./ vals .- 1) .- lookback_time(c, h.z) # in Gyr
end


"""
    convert_units_solar_metallicity(zs::AbstractMatrix{<:Number}, mass::AbstractVector{<:Number})

Converts an `11×n` mass matrix and an `n`-length mass vector to a metallicity vector of length `n`,
in solar metallicities.
"""
function convert_units_solar_metallicity(zs::AbstractMatrix{<:Number}, mass::AbstractVector{<:Number})
    dropdims(sum(zs[2:11, :]; dims=1); dims=1) ./ (mass .- dropdims(sum(zs; dims=1); dims=1)) ./ 0.0134 # in solar metallicities
end
