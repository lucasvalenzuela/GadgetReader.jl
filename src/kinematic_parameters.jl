total_mass(n::Integer, mass::Number) = n * mass
total_mass(mass::AbstractVector{<:Number}) = sum(mass)
total_mass(pos::AbstractMatrix{<:Number}, mass::Number) = total_mass(size(pos, 2), mass)
total_mass(pos::AbstractMatrix{<:Number}, mass::AbstractVector{<:Number}) = total_mass(mass)


@doc raw"""
    angular_momentum(pos, vel, mass[, ellipsoid])
    angular_momentum(p::Particles[, ellipsoid])

Returns the total angular momentum ``J = \sum_i m_i \mathbf{r}_i \times \mathbf{v}_i`` within the
ellipsoid if given. Values are assumed to be in kpc, km/s, and solar masses.

```julia
using Unitful, UnitfulAstro

g::Galaxy # assume an already read-in galaxy with :full units
sph = Sphere(8u"kpc")

J = angular_momentum(g.stars, sph)
```
"""
function angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Union{Number,AbstractVector{<:Number}},
)
    return dropdims(sum(mass' .* pos × vel; dims=2); dims=2)
end

function angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number},
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return angular_momentum(pos[:, mask], vel[:, mask], mass[mask])
end

function angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Number,
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return angular_momentum(pos[:, mask], vel[:, mask], mass)
end

function angular_momentum(p::Particles, ellipsoid::Union{Ellipsoid,Nothing}=nothing;)
    if isnothing(ellipsoid)
        return angular_momentum(p.pos, p.vel, p.mass)
    else
        return angular_momentum(p.pos, p.vel, p.mass, ellipsoid)
    end
end


@doc raw"""
    specific_angular_momentum(pos, vel, mass; Mtot)
    specific_angular_momentum(pos, vel, mass, ellipsoid)
    specific_angular_momentum(p::Particles[, ellipsoid]; Mtot)

Returns the specific angular momentum ``j = \sum_i m_i \mathbf{r}_i \times \mathbf{v}_i / \sum_i m_i``
within the ellipsoid if given. Values are assumed to be in kpc, km/s, and solar masses.
The total mass can be provided via `Mtot` if already known.

```julia
using Unitful, UnitfulAstro

g::Galaxy # assume an already read-in galaxy with :full units
sph = Sphere(8u"kpc")

j = specific_angular_momentum(g.stars, sph)
```
"""
function specific_angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Union{Number,AbstractVector{<:Number}};
    Mtot::Number=total_mass(pos, mass),
)
    return angular_momentum(pos, vel, mass) ./ Mtot
end

function specific_angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Number,
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return specific_angular_momentum(pos[:, mask], vel[:, mask], mass)
end

function specific_angular_momentum(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number},
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return specific_angular_momentum(pos[:, mask], vel[:, mask], mass[mask])
end

function specific_angular_momentum(
    p::Particles,
    ellipsoid::Union{Ellipsoid,Nothing}=nothing;
    Mtot::Union{Number,Nothing}=nothing,
)
    if isnothing(ellipsoid)
        if isnothing(Mtot)
            return specific_angular_momentum(p.pos, p.vel, p.mass)
        else
            return specific_angular_momentum(p.pos, p.vel, p.mass; Mtot)
        end
    else
        return specific_angular_momentum(p.pos, p.vel, p.mass, ellipsoid)
    end
end


@doc raw"""
    b_value(pos, vel, mass[, Mtot])
    b_value(pos, vel, mass, ellipsoid)
    b_value(p::Particles[, ellipsoid])

Returns the b-value ``b = \log(j / (\mathrm{km/s})) - 2/3 \log(M_\mathrm{tot}/\mathrm{M}_\odot)``,
where ``j`` is the [`specific_angular_momentum`](@ref). Values are assumed to be in kpc, km/s, and
solar masses.

```julia
using Unitful, UnitfulAstro

g::Galaxy # assume an already read-in galaxy with :full units
sph = Sphere(8u"kpc")

b = b_value(g.stars, sph)
```
"""
function b_value(
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real},
    mass::Union{Number,AbstractVector{<:Real}},
    Mtot::Real,
)
    return log10(norm(specific_angular_momentum(pos, vel, mass; Mtot))) - 2 // 3 * log10(Mtot)
end

function b_value(
    pos::AbstractMatrix{<:Quantity},
    vel::AbstractMatrix{<:Quantity},
    mass::Union{Number,AbstractVector{<:Quantity}},
    Mtot::Quantity,
)
    return log10(ustrip(u"kpc * km/s", norm(specific_angular_momentum(pos, vel, mass; Mtot)))) -
           2 // 3 * log10(ustrip(u"Msun", Mtot))
end

function b_value(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Union{Number,AbstractVector{<:Number}},
)
    return b_value(pos, vel, mass, total_mass(pos, mass))
end

function b_value(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Number,
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return b_value(pos[:, mask], vel[:, mask], mass)
end

function b_value(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number},
    ellipsoid::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, ellipsoid)
    return b_value(pos[:, mask], vel[:, mask], mass[mask])
end

b_value(p::Particles, ellipsoid::Nothing=nothing) = b_value(p.pos, p.vel, p.mass)
b_value(p::Particles, ellipsoid::Ellipsoid) = b_value(p.pos, p.vel, p.mass, ellipsoid)


@doc raw"""
    circular_velocity(Mtot, r)

Returns the circular velocity ``v_\mathrm{circ} = \sqrt{G M(<r) / r}``. Values are assumed to be in
solar masses and kpc.
"""
function circular_velocity(Mtot::T1, r::T2) where {T1<:AbstractFloat,T2<:AbstractFloat}
    return convert(promote_type(T1, T2), 1e-3sqrt(1.3271244e20Mtot / (3.0856775814913674e19r))) # in km/s
end
circular_velocity(Mtot::Integer, r::T) where {T<:AbstractFloat} = circular_velocity(convert(T, Mtot), r)
circular_velocity(Mtot::T, r::Integer) where {T<:AbstractFloat} = circular_velocity(Mtot, convert(T, r))
circular_velocity(Mtot::Integer, r::Integer) = circular_velocity(float(Mtot), float(r))

function circular_velocity(Mtot::Quantity, r::Quantity)
    vcirc = sqrt(ustrip(u"Msun", Mtot) * u"GMsun" / r)

    # convert to float type of vcirc (conversion to km/s defaults to Float64)
    return convert(typeof(vcirc).parameters[1], ustrip(u"km/s", vcirc)) * u"km/s"
end


@doc raw"""
    spin_parameter(j, Mtot)
    spin_parameter(pos, vel, mass, Mtot, r)
    spin_parameter(pos, vel, mass, Mtot, sphere)
    spin_parameter(p::Particles, Mtot, r; mask=nothing)

Returns the spin parameter ``\lambda(r) = j / (\sqrt{2} r v_\mathrm{circ}(r))``, where ``j`` is
the [`specific_angular_momentum`](@ref). The spin parameter is computed for only one particle component.

A mask can be applied to filter the particles.
"""
function spin_parameter(j::AbstractVector{<:Number}, Mtot::Number, r::Number)
    return norm(j) / (√2 * r * circular_velocity(Mtot, r))
end

function spin_parameter(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Union{Number,AbstractVector{<:Number}},
    Mtot::Number,
    r::Number,
)
    return spin_parameter(specific_angular_momentum(pos, vel, mass), Mtot, r)
end

function spin_parameter(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::Number,
    Mtot::Number,
    sphere::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, sphere)
    return spin_parameter(pos[:, mask], vel[:, mask], mass, Mtot, sphere.radius)
end

function spin_parameter(
    pos::AbstractMatrix{<:Number},
    vel::AbstractMatrix{<:Number},
    mass::AbstractVector{<:Number},
    Mtot::Number,
    sphere::Ellipsoid,
)
    mask = ellipsoidal_mask(pos, sphere)
    return spin_parameter(pos[:, mask], vel[:, mask], mass[mask], Mtot, sphere.radius)
end

function spin_parameter(p::Particles, Mtot::Number, r::Number; mask::Union{AbstractVector,Nothing}=nothing)
    if isnothing(mask)
        return spin_parameter(p.pos, p.vel, p.mass, Mtot, Sphere(r))
    end
    return @views spin_parameter(
        p.pos[:, mask],
        p.vel[:, mask],
        p.mass isa Number ? p.mass : p.mass[mask],
        Mtot,
        Sphere(r),
    )
end

@doc raw"""
    spin_parameter(g::AbstractGalaxy, Mtot, r)

Returns the total spin parameter ``\lambda(r) = j / (\sqrt{2} r v_\mathrm{circ}(r))``, where ``j`` is
the [`specific_angular_momentum`](@ref). The spin parameter is computed for all particle components.
"""
function spin_parameter(g::AbstractGalaxy, Mtot::Number, r::Number)
    sph = Sphere(r)
    J = angular_momentum(g.stars, sph) .+ angular_momentum(g.gas, sph) .+ angular_momentum(g.dm, sph)
    return spin_parameter(J ./ Mtot, Mtot, r)
end
