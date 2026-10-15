"""
    Multibody2D

Minimal **symbolic planar multibody kinematics and energies**, built on
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) (a lightweight replacement for
RigidBodyDynamics.jl when the model must stay symbolic).

# Concepts

* A [`MultibodySystem`](@ref) holds the generalized coordinates `q`, their velocities `q̇`
  (both vectors of symbolic variables) and the gravity vector.
* A [`Frame`](@ref) is a 2D reference frame given by a translation vector and a rotation angle
  with respect to its parent frame (the `world` frame, or another frame).
* A [`Point`](@ref) is a location in a frame, optionally carrying a rigid body
  (mass, center of mass, moment of inertia about the center of mass).

Translations and angles may be numbers, symbolic parameters, or expressions of the generalized
coordinates `q`. **Everything time varying must be written in terms of `q`**: the velocity of any
quantity is obtained as `∂(·)/∂q * q̇`, there is no separate time variable.

# Example (double pendulum)

```julia
using Symbolics
include("Multibody2D.jl")
using .Multibody2D

@variables q[1:2] qd[1:2] l1 l2 lc1 lc2 m1 m2 I1 I2 g
sys = MultibodySystem(collect(q), collect(qd); gravity = g)

link1 = add_frame!(sys, :link1, sys.world, [0, 0], q[1])
link2 = add_frame!(sys, :link2, link1, [l1, 0], q[2])
add_point!(sys, :body1, link1, [lc1, 0]; mass = m1, inertia = I1)
add_point!(sys, :body2, link2, [lc2, 0]; mass = m2, inertia = I2)
tip = add_point!(sys, :tip, link2, [l2, 0])

position(tip)             # tip position in the world frame
velocity(sys, tip)        # tip velocity in the world frame
D = mass_matrix(sys)      # kinetic energy is T = 1/2 q̇ᵀ D(q) q̇
V = potential_energy(sys)

# Dynamics in normal form, ẋ = f(x) + g(x) τ with x = (q, v_a, p_u): only q2 is actuated
# (acrobot), so p_u = p₁ = ∂L/∂q̇₁ replaces q̇₁ and its dynamics do not involve τ.
nf = normal_form(sys, selection_matrix(2, [2]))
nf.state, nf.xdot         # x = [q1, q2, q̇2, p₁] and ẋ(x, τ)
nf.drift, nf.gain         # τ-independent part and input matrix of ẋ
nf.velocities             # q̇ as a function of x
```
"""
module Multibody2D

using LinearAlgebra
using Symbolics
using Symbolics: Num, unwrap, wrap, iscall, operation, arguments

const SU = Symbolics.SymbolicUtils

export Frame, Point, MultibodySystem
export add_frame!, add_point!, relative_transform, rotation_matrix
export position, velocity
export mass_matrix, kinetic_energy, potential_energy, lagrangian
export bias_vector, selection_matrix, normal_form, NormalForm
export trig_reduce

# ------------------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------------------

_num(x) = x isa Num ? x : Num(x)

function _vec2(x, what::AbstractString)
    length(x) == 2 || throw(ArgumentError("$what must have 2 entries, got $(length(x))"))
    return Num[_num(x[1]), _num(x[2])]
end

"""
    rotation_matrix(θ)

The 2×2 rotation matrix `R(θ) = [cos θ -sin θ; sin θ cos θ]`, as a `Matrix{Num}`.
"""
function rotation_matrix(θ)
    θ = _num(θ)
    c, s = cos(θ), sin(θ)
    return Num[c -s; s c]
end

# Transpose of R(θ) written out explicitly, so that no `cos(-θ)` / `sin(-θ)` appears.
function _rotation_matrix_transpose(θ)
    θ = _num(θ)
    c, s = cos(θ), sin(θ)
    return Num[c s; -s c]
end

# ------------------------------------------------------------------------------------------
# Frames
# ------------------------------------------------------------------------------------------

"""
    Frame

A 2D reference frame. Its origin sits at `translation` (expressed in the `parent` frame) and its
axes are rotated by `angle` (counter-clockwise, radians) with respect to the axes of `parent`.
A point with coordinates `x` in this frame has coordinates `translation + R(angle) * x` in the
parent frame. The root of every frame tree is a frame with `parent === nothing`
(`sys.world` for a [`MultibodySystem`](@ref)).

Create frames with [`add_frame!`](@ref).
"""
struct Frame
    name::Symbol
    parent::Union{Nothing, Frame}
    translation::Vector{Num}
    angle::Num
end

function Frame(name::Symbol, parent::Union{Nothing, Frame}, translation, angle)
    return Frame(name, parent, _vec2(translation, "translation"), _num(angle))
end

# The root ("world") frame of a tree.
Frame(name::Symbol = :world) = Frame(name, nothing, Num[0, 0], Num(0))

function Base.show(io::IO, f::Frame)
    if f.parent === nothing
        return print(io, "Frame(:", f.name, ", root)")
    end
    return print(
        io,
        "Frame(:",
        f.name,
        " in :",
        f.parent.name,
        ", translation = ",
        f.translation,
        ", angle = ",
        f.angle,
        ")",
    )
end

function _root(f::Frame)
    while f.parent !== nothing
        f = f.parent
    end
    return f
end

# Frames from `f` up to (and including) the root.
function _ancestors(f::Frame)
    chain = Frame[f]
    while f.parent !== nothing
        f = f.parent
        push!(chain, f)
    end
    return chain
end

# ------------------------------------------------------------------------------------------
# Rigid planar transforms: x_to = t + R(a) * x_from
# ------------------------------------------------------------------------------------------

struct _Transform
    t::Vector{Num}
    a::Num
end

_identity() = _Transform(Num[0, 0], Num(0))

# Apply `T1` first and then `T2`.
function _compose(T2::_Transform, T1::_Transform)
    return _Transform(T2.t + rotation_matrix(T2.a) * T1.t, T2.a + T1.a)
end

function _inverse(T::_Transform)
    return _Transform(-(_rotation_matrix_transpose(T.a) * T.t), -T.a)
end

_apply(T::_Transform, x) = T.t + rotation_matrix(T.a) * x

_transform_to_parent(f::Frame) = _Transform(f.translation, f.angle)

# Transform from the coordinates of `f` to the coordinates of its ancestor `anc`.
function _transform_to_ancestor(f::Frame, anc::Frame)
    T = _identity()
    node = f
    while node !== anc
        node.parent === nothing &&
            error("frame :$(anc.name) is not an ancestor of :$(f.name)")
        T = _compose(_transform_to_parent(node), T)
        node = node.parent
    end
    return T
end

function _transform(from::Frame, to::Frame)
    _root(from) === _root(to) || throw(
        ArgumentError("frames :$(from.name) and :$(to.name) belong to different trees"),
    )
    anc_to = _ancestors(to)
    lca = first(a for a in _ancestors(from) if any(b -> b === a, anc_to))
    up = _transform_to_ancestor(from, lca)
    down = _transform_to_ancestor(to, lca)
    return _compose(_inverse(down), up)
end

"""
    relative_transform(from::Frame, to::Frame) -> (t, a)

The planar transform from the coordinates of `from` to the coordinates of `to`: a point with
coordinates `x` in `from` has coordinates `t + R(a) * x` in `to`. The path goes through the
lowest common ancestor of the two frames, so the angle `a` is expressed with the fewest terms.
"""
function relative_transform(from::Frame, to::Frame)
    T = _transform(from, to)
    return (T.t, T.a)
end

# ------------------------------------------------------------------------------------------
# Points and bodies
# ------------------------------------------------------------------------------------------

# Rigid body rigidly attached to a frame; `com` is in the coordinates of that frame.
struct _Body
    mass::Num
    com::Vector{Num}
    inertia::Num
end

"""
    Point

A point of the system, located at `position` (coordinates in `frame`). If `body` is not
`nothing` the point carries a rigid body of given mass, center of mass `com` (coordinates in
`frame`) and moment of inertia about the center of mass, rigidly attached to `frame`.

Create points with [`add_point!`](@ref).
"""
struct Point
    name::Symbol
    frame::Frame
    position::Vector{Num}
    body::Union{Nothing, _Body}
end

function Base.show(io::IO, p::Point)
    print(io, "Point(:", p.name, " in :", p.frame.name, ", position = ", p.position)
    if p.body !== nothing
        print(
            io,
            ", mass = ",
            p.body.mass,
            ", com = ",
            p.body.com,
            ", inertia = ",
            p.body.inertia,
        )
    end
    return print(io, ")")
end

# ------------------------------------------------------------------------------------------
# System
# ------------------------------------------------------------------------------------------

"""
    MultibodySystem(q, q̇; gravity = 9.81)

A planar multibody system with generalized coordinates `q` and velocities `q̇` (equal-length
vectors of symbolic variables, e.g. `collect(q)` for `@variables q[1:n]`).

`gravity` is either a scalar `g` (gravity acceleration `[0, -g]`, i.e. world `y` pointing up)
or a 2-vector giving the gravity acceleration in the world frame. It may be symbolic.

`sys.world` is the world frame; build the rest of the system with [`add_frame!`](@ref) and
[`add_point!`](@ref), and look entities up by name with `sys[:name]`.
"""
struct MultibodySystem
    q::Vector{Num}
    qd::Vector{Num}
    gravity::Vector{Num}
    world::Frame
    frames::Vector{Frame}   # insertion order, `world` excluded
    points::Vector{Point}   # insertion order
end

function MultibodySystem(q::AbstractVector, qd::AbstractVector; gravity = 9.81)
    length(q) == length(qd) || throw(DimensionMismatch("q and q̇ must have the same length"))
    isempty(q) && throw(ArgumentError("a system needs at least one generalized coordinate"))
    g = gravity isa AbstractVector ? _vec2(gravity, "gravity") : Num[0, -_num(gravity)]
    return MultibodySystem(
        Num[_num(x) for x in q],
        Num[_num(x) for x in qd],
        g,
        Frame(:world),
        Frame[],
        Point[],
    )
end

function Base.show(io::IO, sys::MultibodySystem)
    return print(
        io,
        "MultibodySystem(",
        length(sys.q),
        " coordinates, ",
        length(sys.frames),
        " frames, ",
        length(sys.points),
        " points)",
    )
end

function _name_taken(sys::MultibodySystem, name::Symbol)
    return name === sys.world.name ||
           any(f -> f.name === name, sys.frames) ||
           any(p -> p.name === name, sys.points)
end

function _check_new_name(sys::MultibodySystem, name::Symbol)
    _name_taken(sys, name) && throw(ArgumentError("the name :$name is already used"))
    return nothing
end

function _check_frame(sys::MultibodySystem, f::Frame)
    _root(f) === sys.world ||
        throw(ArgumentError("frame :$(f.name) does not belong to this system"))
    return nothing
end

"""
    Base.getindex(sys::MultibodySystem, name::Symbol)

The frame or point called `name` (`sys[:world]` is the world frame).
"""
function Base.getindex(sys::MultibodySystem, name::Symbol)
    name === sys.world.name && return sys.world
    for f in sys.frames
        f.name === name && return f
    end
    for p in sys.points
        p.name === name && return p
    end
    return throw(KeyError(name))
end

"""
    add_frame!(sys, name, parent, translation, angle) -> Frame

Create a frame called `name`, whose origin is at `translation` (2-vector, in the coordinates of
`parent`) and whose axes are rotated by `angle` with respect to those of `parent`.
`parent` is `sys.world` or a frame previously created in `sys`. Translation entries and the
angle may be numbers, symbolic parameters or expressions of the generalized coordinates.
"""
function add_frame!(sys::MultibodySystem, name::Symbol, parent::Frame, translation, angle)
    _check_new_name(sys, name)
    _check_frame(sys, parent)
    f = Frame(name, parent, translation, angle)
    push!(sys.frames, f)
    return f
end

"""
    add_point!(sys, name, frame, position; mass = nothing, com = position, inertia = 0) -> Point

Create a point called `name` at `position` (2-vector, coordinates in `frame`).

Without `mass` the point is purely kinematic (e.g. a foot or a joint location). With `mass`, a
rigid body is attached to `frame` through this point:

* `mass`: scalar mass of the body;
* `com`: position of the body's center of mass, a 2-vector in the coordinates of `frame`
  (default: the point itself, i.e. a point mass at `position`);
* `inertia`: scalar moment of inertia of the body about its center of mass (default `0`).

All quantities may be numbers or symbolic.
"""
function add_point!(
    sys::MultibodySystem,
    name::Symbol,
    frame::Frame,
    position;
    mass = nothing,
    com = nothing,
    inertia = nothing,
)
    _check_new_name(sys, name)
    _check_frame(sys, frame)
    pos = _vec2(position, "position")
    if mass === nothing
        (com === nothing && inertia === nothing) ||
            throw(ArgumentError("`com` and `inertia` require `mass`"))
        body = nothing
    else
        body = _Body(
            _num(mass),
            com === nothing ? copy(pos) : _vec2(com, "com"),
            inertia === nothing ? Num(0) : _num(inertia),
        )
    end
    p = Point(name, frame, pos, body)
    push!(sys.points, p)
    return p
end

# ------------------------------------------------------------------------------------------
# Trigonometric simplification
# ------------------------------------------------------------------------------------------

# Symbolics' `simplify` knows `sin² + cos² = 1` but not product-to-sum identities, so products such as
# `cos(q1) cos(q1 + q2) + sin(q1) sin(q1 + q2)` (which Jacobian products produce all the time) are
# left alone. `trig_reduce` rewrites a polynomial in `sin`/`cos` as a *linear* combination of
# `sin(arg)`/`cos(arg)` with canonical arguments, where like terms (and hence all Pythagorean
# identities) cancel exactly.

function _constant(x)
    x isa Number && return x
    v = Symbolics.value(x)
    v isa Number && return v
    if isdefined(SU, :isconst) && SU.isconst(v)
        return SU.unwrap_const(v)
    end
    return nothing
end

_is_trig(e) = iscall(e) && (operation(e) === sin || operation(e) === cos)

# Split `c * rest` into the numeric `c` and `rest` (`nothing` for a pure number).
function _numeric_coefficient(t)
    c = _constant(t)
    c === nothing || return (c, nothing)
    if iscall(t) && operation(t) === (*)
        args = arguments(t)
        c = _constant(args[1])
        if c !== nothing
            rest = length(args) == 2 ? args[2] : prod(args[2:end])
            return (c, rest)
        end
    end
    return (1, t)
end

# Canonical form of the argument `a` of a trig function, and the sign change it induced:
# `cos(-a) = cos(a)` and `sin(-a) = -sin(a)`, so `a` and `-a` must share one representative.
function _canonical_angle(a)
    a = unwrap(expand(_num(a)))
    terms = iscall(a) && operation(a) === (+) ? arguments(a) : [a]
    best_key, best_coef = nothing, 1
    for t in terms
        c, rest = _numeric_coefficient(t)
        rest === nothing && continue
        key = string(rest)
        if best_key === nothing || key < best_key
            best_key, best_coef = key, c
        end
    end
    flip = best_key !== nothing && best_coef < 0
    return flip ? (unwrap(expand(-wrap(a))), true) : (a, false)
end

const _ONE_KEY = (:one, nothing)

function _add_term!(acc, kind::Symbol, arg, coef)
    arg_c, flipped = _canonical_angle(arg)
    zero_arg = _constant(arg_c) !== nothing && iszero(_constant(arg_c))
    if zero_arg
        kind === :sin && return nothing
        key = _ONE_KEY
    else
        key = (kind, arg_c)
        if kind === :sin && flipped
            coef = -coef
        end
    end
    acc[key] = get(acc, key, Num(0)) + coef
    return nothing
end

function _product_to_sum(k1::Symbol, a1, k2::Symbol, a2)
    s = a1 + a2
    d = a1 - a2
    half = Num(1 // 2)
    if k1 === :cos && k2 === :cos
        return ((:cos, d, half), (:cos, s, half))
    elseif k1 === :sin && k2 === :sin
        return ((:cos, d, half), (:cos, s, -half))
    elseif k1 === :sin && k2 === :cos
        return ((:sin, s, half), (:sin, d, half))
    else
        return ((:sin, s, half), (:sin, d, -half))
    end
end

# Collects `scale * e` into `acc`, a dictionary (kind, arg) => coefficient.
function _collect!(acc, e, scale)
    if iscall(e) && operation(e) === +
        for a in arguments(e)
            _collect!(acc, a, scale)
        end
        return nothing
    end
    coef = Num(1)
    trigs = Tuple{Symbol, Any}[]
    factors = iscall(e) && operation(e) === (*) ? arguments(e) : [e]
    for f in factors
        if _is_trig(f)
            push!(trigs, (nameof(operation(f)), wrap(only(arguments(f)))))
        elseif iscall(f) &&
               operation(f) === (^) &&
               _is_trig(arguments(f)[1]) &&
               (n = _constant(arguments(f)[2])) isa Integer &&
               n > 0
            base = arguments(f)[1]
            for _ in 1:n
                push!(trigs, (nameof(operation(base)), wrap(only(arguments(base)))))
            end
        else
            coef *= _reduce_atom(f)
        end
    end
    # Multiply the trig factors one by one, linearizing as we go.
    terms = Pair{Tuple{Symbol, Any}, Num}[_ONE_KEY => scale * coef]
    for (k2, a2) in trigs
        new_terms = Pair{Tuple{Symbol, Any}, Num}[]
        for ((k1, a1), c) in terms
            if k1 === :one
                push!(new_terms, (k2, a2) => c)
            else
                for (k, a, w) in _product_to_sum(k1, wrap(a1), k2, a2)
                    push!(new_terms, (k, a) => c * w)
                end
            end
        end
        terms = new_terms
    end
    for ((k, a), c) in terms
        if k === :one
            acc[_ONE_KEY] = get(acc, _ONE_KEY, Num(0)) + c
        else
            _add_term!(acc, k, a, c)
        end
    end
    return nothing
end

# A factor that is neither a sum, a product nor a trig function: constants, symbols, and generic
# function calls (whose arguments are reduced recursively).
function _reduce_atom(f)
    # `q[i]` (an element of an array variable) is an atom, its arguments are not expressions.
    if iscall(f) &&
       operation(f) !== (+) &&
       operation(f) !== (*) &&
       operation(f) !== getindex
        try
            args = arguments(f)
            reduced = map(a -> unwrap(trig_reduce(wrap(a))), args)
            if !all(i -> isequal(reduced[i], args[i]), eachindex(args))
                return wrap(operation(f)(reduced...))
            end
        catch
            # Unknown function or exotic argument: keep the atom as it is.
        end
    end
    return _num(wrap(f))
end

"""
    trig_reduce(x)

Rewrite the symbolic expression `x` (or each entry of an array) as a linear combination of
`sin(arg)`/`cos(arg)` with canonical arguments, using the product-to-sum identities.
This is what turns `cos(q1) cos(q1 + q2) + sin(q1) sin(q1 + q2)` into `cos(q2)` and makes
every `sin² + cos²` cancel, which `Symbolics.simplify` alone does not do. Non-trigonometric
structure (parameters, polynomials in `q̇`, ...) is only expanded.
"""
function trig_reduce(x::Num)
    e = unwrap(expand(x))
    acc = Dict{Tuple{Symbol, Any}, Num}()
    _collect!(acc, e, Num(1))
    result = Num(0)
    for key in sort!(collect(keys(acc)); by = k -> (string(k[1]), string(k[2])))
        coef = expand(acc[key])
        iszero(coef) && continue
        kind, arg = key
        trig = kind === :one ? Num(1) : (kind === :sin ? sin(wrap(arg)) : cos(wrap(arg)))
        result += coef * trig
    end
    return result
end

trig_reduce(x::Real) = x
trig_reduce(A::AbstractArray) = map(trig_reduce, A)

_finish(x, simplify::Bool) = simplify ? trig_reduce(x) : x

# ------------------------------------------------------------------------------------------
# Position and velocity
# ------------------------------------------------------------------------------------------

"""
    position(p::Point, frame::Frame = root frame of p) -> Vector{Num}

Coordinates of the point `p` in `frame`. `simplify = true` (default) applies [`trig_reduce`](@ref).
"""
function Base.position(p::Point, frame::Frame = _root(p.frame); simplify::Bool = true)
    t, a = relative_transform(p.frame, frame)
    return _finish.(t + rotation_matrix(a) * p.position, simplify)
end

"""
    velocity(sys, p::Point; relative_to = sys.world, expressed_in = relative_to,
             simplify = true) -> Vector{Num}

Velocity of the point `p`, as seen by an observer attached to the frame `relative_to` (the world
frame by default), with components in the axes of the frame `expressed_in`.

The two keywords are different things: `relative_to` selects *who observes* (the velocity of a
point at rest in `relative_to` is zero), `expressed_in` only selects *which axes* the
components are written in. For instance `velocity(sys, p; expressed_in = link)` is the world
velocity of `p` written in the axes of `link`.

The velocity is `∂(position)/∂q * q̇`, so any time dependence must go through `q`.
"""
function velocity(
    sys::MultibodySystem,
    p::Point;
    relative_to::Frame = sys.world,
    expressed_in::Frame = relative_to,
    simplify::Bool = true,
)
    _check_frame(sys, p.frame)
    _check_frame(sys, relative_to)
    _check_frame(sys, expressed_in)
    x = position(p, relative_to; simplify = false)
    v = Symbolics.jacobian(x, sys.q) * sys.qd
    if expressed_in !== relative_to
        _, a = relative_transform(relative_to, expressed_in)
        v = rotation_matrix(a) * v
    end
    return _finish.(v, simplify)
end

# ------------------------------------------------------------------------------------------
# Energies
# ------------------------------------------------------------------------------------------

# World angle of `frame`, world position of the point `x` given in `frame`.
function _world_pose(frame::Frame, x)
    T = _transform_to_ancestor(frame, _root(frame))
    return (_apply(T, x), T.a)
end

"""
    mass_matrix(sys; simplify = true) -> Matrix{Num}

The symmetric matrix `D(q)` such that the kinetic energy of all bodies of `sys` is
`T = 1/2 q̇ᵀ D(q) q̇`. Each body contributes `m JᵀJ + I J_θᵀ J_θ`, with `J = ∂(center of mass)/∂q`
and `J_θ = ∂(frame angle)/∂q`, both taken in the world frame.

`simplify = true` (default) reduces every entry with [`trig_reduce`](@ref).
"""
function mass_matrix(sys::MultibodySystem; simplify::Bool = true)
    n = length(sys.q)
    D = zeros(Num, n, n)
    for p in sys.points
        p.body === nothing && continue
        b = p.body
        pc, θ = _world_pose(p.frame, b.com)
        J = Symbolics.jacobian(pc, sys.q)
        Jθ = Symbolics.jacobian([θ], sys.q)
        D += b.mass * (J' * J) + b.inertia * (Jθ' * Jθ)
    end
    # Only the upper triangle is simplified; D is symmetric by construction.
    for i in 1:n, j in i:n
        D[i, j] = _finish(expand(D[i, j]), simplify)
        D[j, i] = D[i, j]
    end
    return D
end

"""
    kinetic_energy(sys; simplify = true) -> Num

Kinetic energy `1/2 q̇ᵀ D(q) q̇`, with `D` from [`mass_matrix`](@ref).
"""
function kinetic_energy(sys::MultibodySystem; simplify::Bool = true)
    D = mass_matrix(sys; simplify)
    n = length(sys.qd)
    # Explicit double sum: `q̇' * D * q̇` would go through Symbolics' array machinery.
    T = sum(sys.qd[i] * D[i, j] * sys.qd[j] for i in 1:n, j in 1:n) / 2
    return _finish(expand(T), simplify)
end

"""
    potential_energy(sys; simplify = true) -> Num

Gravitational potential energy `-Σ m gᵀ r_com` of all bodies of `sys`, with `g` the gravity
acceleration of the system and `r_com` the world position of each center of mass.
"""
function potential_energy(sys::MultibodySystem; simplify::Bool = true)
    V = Num(0)
    for p in sys.points
        p.body === nothing && continue
        pc, _ = _world_pose(p.frame, p.body.com)
        V -= p.body.mass * sum(sys.gravity .* pc)
    end
    return _finish(expand(V), simplify)
end

"""
    lagrangian(sys; simplify = true) -> Num

The Lagrangian `L = T - V` (kinetic minus potential energy) of `sys`.
"""
function lagrangian(sys::MultibodySystem; simplify::Bool = true)
    return _finish(
        expand(kinetic_energy(sys; simplify) - potential_energy(sys; simplify)),
        simplify,
    )
end

# ------------------------------------------------------------------------------------------
# Equations of motion
# ------------------------------------------------------------------------------------------

# Returns (D, h, f) with `f = ∂L/∂q` and `h = ∂L/∂q - Ḋ q̇` (all functions of q and q̇), so that
# the Euler-Lagrange equations `d/dt(D q̇) - ∂L/∂q = Bτ` read `D q̈ = h + Bτ`.
function _lagrangian_terms(sys::MultibodySystem; simplify::Bool)
    n = length(sys.q)
    q, qd = sys.q, sys.qd
    D = mass_matrix(sys; simplify)
    V = potential_energy(sys; simplify)
    dD = [map(x -> Symbolics.derivative(x, q[k]), D) for k in 1:n]
    f = Vector{Num}(undef, n)
    h = Vector{Num}(undef, n)
    for i in 1:n
        # ∂L/∂q_i = 1/2 q̇ᵀ (∂D/∂q_i) q̇ - ∂V/∂q_i
        quad = sum(dD[i][j, k] * qd[j] * qd[k] for j in 1:n, k in 1:n) / 2
        # (Ḋ q̇)_i = Σ_jk (∂D_ij/∂q_k) q̇_k q̇_j
        Ddot = sum(dD[k][i, j] * qd[k] * qd[j] for k in 1:n, j in 1:n)
        f[i] = _finish(expand(quad - Symbolics.derivative(V, q[i])), simplify)
        h[i] = _finish(expand(quad - Symbolics.derivative(V, q[i]) - Ddot), simplify)
    end
    return D, h, f
end

"""
    bias_vector(sys; simplify = true) -> Vector{Num}

The vector `h(q, q̇) = ∂L/∂q - Ḋ(q) q̇` (Coriolis, centrifugal and gravity terms, sign
included) for which the Euler-Lagrange equations of `sys` read `D(q) q̈ = h(q, q̇) + B τ`.
"""
function bias_vector(sys::MultibodySystem; simplify::Bool = true)
    return _lagrangian_terms(sys; simplify)[2]
end

"""
    selection_matrix(n, actuated) -> Matrix{Num}

The `n × length(actuated)` input matrix `B` in which input `j` acts directly on the generalized
coordinate `actuated[j]` (column `j` is the `actuated[j]`-th Euclidean basis vector). The rows of
the unactuated coordinates are zero.
"""
function selection_matrix(n::Integer, actuated::AbstractVector{<:Integer})
    allunique(actuated) || throw(ArgumentError("`actuated` has repeated indices"))
    all(i -> 1 <= i <= n, actuated) ||
        throw(ArgumentError("`actuated` indices must lie in 1:$n"))
    B = zeros(Num, n, length(actuated))
    for (j, i) in enumerate(actuated)
        B[i, j] = 1
    end
    return B
end

# Solves `A X = b` by Gaussian elimination *without pivoting*, which is legitimate (and keeps the
# symbolic expressions free of case distinctions) because `A` is always symmetric positive
# definite here: a block of the mass matrix or a Schur complement of it.
function _solve(A::AbstractMatrix, b::AbstractVecOrMat)
    n = size(A, 1)
    n == 0 && return copy(b)
    M = Num[A[i, j] for i in 1:n, j in 1:n]
    X = Num[x for x in reshape(b, n, :)]
    for k in 1:n, i in (k + 1):n
        r = M[i, k] / M[k, k]
        for j in (k + 1):n
            M[i, j] -= r * M[k, j]
        end
        X[i, :] .-= r .* X[k, :]
    end
    for k in n:-1:1
        for j in (k + 1):n
            X[k, :] .-= M[k, j] .* X[j, :]
        end
        X[k, :] ./= M[k, k]
    end
    return b isa AbstractVector ? vec(X) : X
end

# Schur complement `S = D_aa - D_au D_uu⁻¹ D_ua`, returned as `(N, den)` with `S = N / den`. For
# one or two unactuated coordinates the adjugate/determinant of `D_uu` keeps `N` a polynomial in
# trigonometric functions (which `trig_reduce` compacts) and `S⁻¹ r = den N⁻¹ r` has no nested
# fractions. Larger blocks fall back to Gaussian elimination (`den = 1`).
function _schur(D_aa, D_au, D_ua, D_uu, simplify::Bool)
    nu = size(D_uu, 1)
    if nu == 1 || nu == 2
        adj =
            nu == 1 ? fill(Num(1), 1, 1) :
            Num[D_uu[2, 2] -D_uu[1, 2]; -D_uu[2, 1] D_uu[1, 1]]
        den = nu == 1 ? D_uu[1, 1] : D_uu[1, 1] * D_uu[2, 2] - D_uu[1, 2] * D_uu[2, 1]
        N = D_aa * den - D_au * adj * D_ua
        return _finish.(expand.(N), simplify), _finish(expand(den), simplify)
    end
    return D_aa - D_au * _solve(D_uu, D_ua), Num(1)
end

"""
    NormalForm

Control-affine state-space model `ẋ = drift + gain * τ` of a [`MultibodySystem`](@ref), with state
`x = (q, v_a, p_u)` (all coordinates, the actuated velocities and the unactuated generalized
momenta). Built by [`normal_form`](@ref). Fields:

* `state::Vector{Num}`: `x = [q; v_a; p_u]`;
* `input::Vector{Num}`: the inputs `τ`;
* `actuated`, `unactuated`: indices of the actuated / unactuated generalized coordinates;
* `velocities::Vector{Num}`: the full `q̇` as a function of the state (`q̇_a = v_a`,
  `q̇_u = D_uu⁻¹ (p_u - D_ua v_a)`);
* `xdot::Vector{Num}`: `ẋ = drift + gain * input`;
* `drift::Vector{Num}`, `gain::Matrix{Num}`: the `τ`-independent part and the input matrix of `ẋ`.
"""
struct NormalForm
    state::Vector{Num}
    input::Vector{Num}
    actuated::Vector{Int}
    unactuated::Vector{Int}
    velocities::Vector{Num}
    xdot::Vector{Num}
    drift::Vector{Num}
    gain::Matrix{Num}
end

function Base.show(io::IO, nf::NormalForm)
    return print(
        io,
        "NormalForm(",
        length(nf.state),
        " states, ",
        length(nf.input),
        " inputs; actuated = ",
        nf.actuated,
        ", unactuated = ",
        nf.unactuated,
        ")",
    )
end

"""
    normal_form(sys, B; input = nothing, momenta = nothing, simplify = true) -> NormalForm

Dynamics of `sys` in normal form. `B` is the `n × m` input matrix of the Euler-Lagrange equations

    d/dt(∂L/∂q̇_i) - ∂L/∂q_i = B_i τ,        L = T - V,

typically built with [`selection_matrix`](@ref). A coordinate whose row of `B` is zero is
*unactuated*; every other coordinate is *actuated*. Writing `a` / `u` for the two index sets, the
state is `x = (q, v_a, p_u)` with `v_a = q̇_a` and `p_u = (∂L/∂q̇)_u = (D q̇)_u`, and

    q̇_a = v_a,        q̇_u = D_uu⁻¹ (p_u - D_ua v_a),
    ṗ_u = ∂L/∂q_u,
    S v̇_a = h_a - D_au D_uu⁻¹ h_u + B_a τ,        S = D_aa - D_au D_uu⁻¹ D_ua,

with `h = ∂L/∂q - Ḋ q̇` (see [`bias_vector`](@ref)) evaluated at the `q̇` above. Since `B_u = 0`,
the momentum dynamics `ṗ_u` do not depend on `τ`.

The velocities of the actuated coordinates are the entries `sys.qd[actuated]`. `input` and
`momenta` optionally supply the symbolic variables for `τ` and `p_u`; by default `τ₁, τ₂, …` and
`p_i` (`i` the coordinate index) are created. `simplify = true` applies [`trig_reduce`](@ref) to
`D` and `h`; the final fractions (through `D_uu⁻¹` and `S⁻¹`) are left as built.
"""
function normal_form(
    sys::MultibodySystem,
    B::AbstractMatrix;
    input = nothing,
    momenta = nothing,
    simplify::Bool = true,
)
    n = length(sys.q)
    size(B, 1) == n || throw(
        DimensionMismatch(
            "B must have one row per generalized coordinate ($n), got $(size(B, 1))",
        ),
    )
    m = size(B, 2)
    Bm = Num[_num(B[i, j]) for i in 1:n, j in 1:m]
    act = [i for i in 1:n if !all(iszero, Bm[i, :])]
    un = [i for i in 1:n if i ∉ act]

    τ =
        input === nothing ? Num[Symbolics.variable(:τ, j) for j in 1:m] :
        Num[_num(x) for x in input]
    length(τ) == m || throw(DimensionMismatch("`input` must have $m entries"))
    p_u =
        momenta === nothing ? Num[Symbolics.variable(:p, i) for i in un] :
        Num[_num(x) for x in momenta]
    length(p_u) == length(un) || throw(
        DimensionMismatch(
            "`momenta` must have $(length(un)) entries (one per unactuated coordinate)",
        ),
    )

    D, h, f = _lagrangian_terms(sys; simplify)
    v_a = sys.qd[act]
    D_aa, D_au, D_ua, D_uu = D[act, act], D[act, un], D[un, act], D[un, un]

    # q̇ as a function of the state, and the substitution that realizes it.
    qdot = copy(sys.qd)
    qdot[un] = _solve(D_uu, p_u - D_ua * v_a)
    subs = Dict(sys.qd[i] => qdot[i] for i in un)
    on_state(x) = Symbolics.substitute(x, subs)

    pdot_u = on_state.(f[un])
    h_a, h_u = on_state.(h[act]), on_state.(h[un])
    N, den = _schur(D_aa, D_au, D_ua, D_uu, simplify)   # S = N / den
    vdot_drift = _solve(N, den * (h_a - D_au * _solve(D_uu, h_u)))
    vdot_gain = _solve(N, den * Bm[act, :])

    state = [sys.q; v_a; p_u]
    drift = [qdot; vdot_drift; pdot_u]
    gain = [zeros(Num, n, m); vdot_gain; zeros(Num, length(un), m)]
    return NormalForm(state, τ, act, un, qdot, drift + gain * τ, drift, gain)
end

end # module
