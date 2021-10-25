"""
    convert_units!(p::Particles, h::AbstractGadgetHeader, units::Symbol=:full)

Converts the units of the particles properties (see [`Particles`](@ref)) from a snapshot or
subfind header in-place.
Use `:full` for units as `Unitful` quantities, `:physical` for values converted
to physical units (kpc, km/s, solar metallicities etc.), or `:sim` for values in simulation units.
The individual metallicities in "Zs" are converted to a one-dimensional `Vector` of values in solar
metallicities in any case.
"""
function convert_units!(p::Particles, h::AbstractGadgetHeader, units::Symbol=:full)
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

function convert_to_unitful!(
    a::AbstractArray{T},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
) where {T<:Real}
    haskey(simunits, prop) || return a

    factor_func, u = simunits[prop]

    if !isnothing(factor_func)
        a .*= factor_func(h)
    end

    if !isnothing(u)
        Q = Quantity{T,dimension(u),typeof(u)}
        return reinterpret(Q, a)
    end

    return a
end

function convert_to_physical!(
    a::AbstractArray{<:Real},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return a

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        a .*= factor_func(h)
    end

    return a
end

function convert_to_sim!(
    a::AbstractArray{<:Real},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return a

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        a ./= factor_func(h)
    end

    return a
end

function convert_to_sim!(
    a::AbstractArray{Quantity},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    astrip = ustrip(a)
    haskey(simunits, prop) || return astrip

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        astrip ./= factor_func(h)
    end

    return astrip
end


function convert_to_unitful(
    a::AbstractArray{T},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
) where {T<:Real}
    haskey(simunits, prop) || return copy(a)

    factor_func, u = simunits[prop]

    if isnothing(factor_func)
        if isnothing(u)
            return copy(a)
        else
            return a .* u
        end
    else
        if isnothing(u)
            return a .* factor_func(h)
        else
            return a .* factor_func(h) .* u
        end
    end
end

function convert_to_physical(
    a::AbstractArray{<:Real},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return copy(a)

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        return a .* factor_func(h)
    end

    return copy(a)
end

function convert_to_sim(
    a::AbstractArray{<:Real},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return copy(a)

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        return a ./ factor_func(h)
    end

    return copy(a)
end

function convert_to_sim(
    a::AbstractArray{Quantity},
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    astrip = ustrip(a)
    haskey(simunits, prop) || return copy(astrip)

    factor_func, _ = simunits[prop]

    if !isnothing(factor_func)
        return astrip ./ factor_func(h)
    end

    return copy(astrip)
end

function convert_to_unitful(a::Real, prop::AbstractString, h::AbstractGadgetHeader, simunits::Dict)
    haskey(simunits, prop) || return a

    factor_func, u = simunits[prop]

    if !isnothing(factor_func)
        a *= factor_func(h)
    end

    if isnothing(u)
        return a
    else
        return a * u
    end
end

function convert_to_physical(
    a::Real,
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return a

    factor_func, _ = simunits[prop]

    if isnothing(factor_func)
        return a
    else
        return a * factor_func(h)
    end
end

function convert_to_sim(
    a::Real,
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    haskey(simunits, prop) || return a

    factor_func, _ = simunits[prop]

    if isnothing(factor_func)
        return a
    else
        return a / factor_func(h)
    end
end

function convert_to_sim(
    a::Quantity,
    prop::AbstractString,
    h::AbstractGadgetHeader,
    simunits::Dict,
)
    astrip = ustrip(a)
    haskey(simunits, prop) || return astrip

    factor_func, _ = simunits[prop]

    if isnothing(factor_func)
        return astrip
    else
        return astrip / factor_func(h)
    end
end



function convert_units_subfind_prop(
    val::Union{Real,AbstractArray{<:Real}},
    prop::AbstractString,
    h::AbstractGadgetHeader,
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


@doc raw"""
    convert_units_age(val::Union{Real,AbstractArray{<:Real}}, h::AbstractGadgetHeader, units::Symbol=:full)
    convert_units_physical_age(val::Real, h)
    convert_units_physical_age(vals::AbstractArray{<:Real}, h)
    convert_units_physical_age!(vals::AbstractArray{<:Real}, h)
    convert_units_full_age(val::Real, h)
    convert_units_full_age(vals::AbstractArray{<:Real}, h)

Converts ages from the scale factor ``a`` via the lookback time according to the cosmology
defined by the snapshot header.
Full returns values in `Unitful` quantities, whereas physical returns the physical value in Gyr.
"""
function convert_units_age(
    val::Union{Real,AbstractArray{<:Real}},
    h::AbstractGadgetHeader,
    units::Symbol=:full,
)
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
