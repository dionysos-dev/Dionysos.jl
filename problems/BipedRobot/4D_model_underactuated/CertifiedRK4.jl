"""
    CertifiedRK4

Certified successors of the **Float64 RK4 map** of a symbolic system, for the grid abstraction of
Dionysos. Nothing here belongs to Dionysos itself: it plugs into the interfaces that Dionysos
already has (a `DiscreteTimeSystemOverApproximation`, the `USER_DEFINED` approximation mode of the
`UniformGridAbstraction` optimizer, and a `get_states_from_set_strict` method for its own image
type), and is built on `CoreMath`, `IntervalArithmetic` 1.x and `Symbolics`.

[`CertifiedRK4GrowthBound`](@ref) builds the Float64 RK4 system and encloses its one-step images.
The system returned by `ST.get_system` is the concrete system to use for simulation and synthesis, on
a nonperiodic Float64 `GridFree`. For every state that the grid quantizer assigns to a source cell,
the successor computed by that system lies in the image that [`overapproximation_map`](@ref) returns
for the cell: the RK4 evolution of the cell, enclosed with interval arithmetic and the roundoff of
the Float64 arithmetic, with the cells that image can quantize into found with the same quantizer.
If an image cannot be enclosed inside the mapped domain, the action is disabled for that cell.
The certificate covers the sampled map: not the continuous flow between samples, nor the truncation
error of RK4 with respect to the ODE.
"""
module CertifiedRK4

import Dionysos: Mapping as MP, System as ST
import CoreMath
import IntervalArithmetic as IA
import LazySets
import MathematicalSystems as MS
import StaticArrays: SVector
import Symbolics

export CertifiedRK4GrowthBound, overapproximation_map

# ------------------------------------------------------------------------------------------
# Deterministic fingerprints (to identify stored data; `Base.hash` is not made for that)
# ------------------------------------------------------------------------------------------

# FNV-1a on 64 bits: the same in every Julia session, version and machine, so a fingerprint can be
# written to a file on one machine and compared on another.
const FINGERPRINT_SEED = 0xcbf29ce484222325
const _FNV_PRIME = 0x00000100000001b3

_fingerprint_step(h::UInt64, b::UInt8) = (h ⊻ b) * _FNV_PRIME

# The 8 bytes of `u`, least significant first, so that the result does not depend on the endianness.
function _fingerprint_word(h::UInt64, u::UInt64)
    for k in 0:7
        h = _fingerprint_step(h, UInt8((u >> (8k)) & 0xff))
    end
    return h
end

"""
    fingerprint(x, h = FINGERPRINT_SEED) -> UInt64

Deterministic 64-bit FNV-1a fingerprint of strings, integers, floats and tuples/arrays of them,
folded into `h` (so that `fingerprint(y, fingerprint(x))` fingerprints `x` then `y`). Every value
carries a type tag and a length, so that `("ab", "c")` and `("a", "bc")` differ. It detects
accidental changes of the data it is applied to (another grid, other dynamics): it is not a
cryptographic hash.
"""
fingerprint(x, h::UInt64 = FINGERPRINT_SEED) = _fingerprint(x, h)

function _fingerprint(s::AbstractString, h::UInt64)
    h = _fingerprint_word(_fingerprint_step(h, 0x04), UInt64(ncodeunits(s)))
    for b in codeunits(s)
        h = _fingerprint_step(h, b)
    end
    return h
end

_fingerprint(x::Integer, h::UInt64) =
    _fingerprint_word(_fingerprint_step(h, 0x01), reinterpret(UInt64, Int64(x)))

_fingerprint(x::AbstractFloat, h::UInt64) =
    _fingerprint_word(_fingerprint_step(h, 0x02), reinterpret(UInt64, Float64(x)))

function _fingerprint(x::Union{Tuple, AbstractArray}, h::UInt64)
    h = _fingerprint_word(_fingerprint_step(h, 0x03), UInt64(length(x)))
    for e in x
        h = _fingerprint(e, h)
    end
    return h
end

"""
    fingerprint_string(h::UInt64) -> String

The 16 hexadecimal digits of a [`fingerprint`](@ref), to store it as text.
"""
fingerprint_string(h::UInt64) = string(h; base = 16, pad = 16)

"""
    fingerprint_file(path) -> String

[`fingerprint_string`](@ref) of the bytes of a file, read in blocks.
"""
function fingerprint_file(path::AbstractString)
    h = FINGERPRINT_SEED
    buffer = Vector{UInt8}(undef, 1 << 20)
    open(path, "r") do io
        while !eof(io)
            n = readbytes!(io, buffer)
            for i in 1:n
                h = _fingerprint_step(h, buffer[i])
            end
        end
    end
    return fingerprint_string(h)
end

# ------------------------------------------------------------------------------------------
# The approximation
# ------------------------------------------------------------------------------------------

"""
    CertifiedRK4GrowthBound(f, x, u, X, U, tstep; nsystem = ST.DEFAULT_NUM_SUBSTEPS)
    CertifiedRK4GrowthBound(system::ConstrainedBlackBoxControlContinuousSystem, tstep; kwargs...)

Conservative successor enclosure of the Float64 RK4 map returned by `ST.get_system`. `f` is branch-free
symbolic dynamics (or a continuous system that can be traced symbolically), `x` and `u` the state and
input variables, `X` and `U` their sets. It is a `DiscreteTimeSystemOverApproximation`, so
`ST.reach_set(approx, cell, u, nothing)` is the enclosure of the image of a box (`nothing` when none
can be certified); use [`overapproximation_map`](@ref) to plug it into the grid abstraction.
"""
struct CertifiedRK4GrowthBound{S, F} <: ST.DiscreteTimeSystemOverApproximation
    system::S
    f::F
    tstep::Float64
    nsystem::Int
    dynamics::String
end

ST.get_system(approx::CertifiedRK4GrowthBound) = approx.system

"""
    KERNEL_VERSION

Version of the certified successor kernel. It is part of [`get_provenance`](@ref): bump it whenever a
change of the kernel can change the transitions it returns, so that cache files and chunks computed
by the older kernel are refused.
"""
const KERNEL_VERSION = 1

"""
    get_provenance(approx::CertifiedRK4GrowthBound) -> NamedTuple

Plain data that identifies **how** the successors of `approx` are computed: the kernel version, the
sampling period, the number of RK4 substeps and a fingerprint of the printed symbolic dynamics.
`BipedAbstraction` stores it with cached transitions and chunks, and refuses transitions that
do not carry the same one.
"""
get_provenance(approx::CertifiedRK4GrowthBound) = (;
    approximation = "CertifiedRK4GrowthBound",
    kernel_version = KERNEL_VERSION,
    tstep = approx.tstep,
    nsystem = approx.nsystem,
    dynamics = approx.dynamics,
)

_to_interval(v::IA.Interval{Float64}) = v
_to_interval(v::Real) = IA.interval(Float64, v)
_point_box(u) = [_to_interval(u[i]) for i in 1:length(u)]
_valid_interval_enclosure(v) = IA.decoration(v) >= IA.dac && IA.isbounded(v)

struct InvalidRK4Evaluation <: Exception end

_checked_primitive(v) = v
_checked_primitive(v::Float64) = isfinite(v) ? v : throw(InvalidRK4Evaluation())
_checked_primitive(v::IA.Interval{Float64}) =
    _valid_interval_enclosure(v) ? v : throw(InvalidRK4Evaluation())

_rk4_operation(op, xs...) = _checked_primitive(op(xs...))
_rk4_shift(x, k, scale) = _rk4_operation.(+, x, _rk4_operation.(*, k, scale))

# Correct rounding keeps the Float64 point result inside an outward-rounded
# interval for the same real operation. Base transcendental functions do not
# offer this guarantee uniformly.
struct RoundedElementary{F}
    f::F
end

# If each scalar result is rounded to nearest, it remains inside its
# outward-rounded interval. Keeping the same operation graph for floats and
# intervals applies this fact at every stage, including roundoff in the slopes.
function _rk4_stages(f, x, u, h)
    k1 = f(x, u)
    y2 = _rk4_shift(x, k1, h / 2)
    k2 = f(y2, u)
    y3 = _rk4_shift(x, k2, h / 2)
    k3 = f(y3, u)
    y4 = _rk4_shift(x, k3, h)
    k4 = f(y4, u)
    sum12 = _rk4_shift(k1, k2, 2.0)
    sum123 = _rk4_shift(sum12, k3, 2.0)
    endpoint = _rk4_shift(x, _rk4_operation.(+, sum123, k4), h / 6)
    return (x, y2, y3, y4), (k1, k2, k3, k4), endpoint
end

function _try_rk4_stages(f, x, u, h)
    try
        return _rk4_stages(f, x, u, h)
    catch err
        err isa InvalidRK4Evaluation || rethrow()
        return nothing
    end
end

function _roundoff_radius(value)
    # Bound the rounding of a real result anywhere in this interval by the
    # largest adjacent Float64 gap at its endpoints.
    _valid_interval_enclosure(value) || return IA.entireinterval(Float64)
    lo, hi = IA.inf(value), IA.sup(value)
    gap = max(lo - prevfloat(lo), nextfloat(hi) - hi)
    return isfinite(gap) ? _to_interval(gap) : IA.entireinterval(Float64)
end

# Propagate the error between a rounded program and its real operation graph.
# Each generated arithmetic call has a point and an interval-error method.
struct RK4ErrorValue <: Real
    range::IA.Interval{Float64}
    error::Float64
end

_error_value(x::RK4ErrorValue) = x
_error_value(x::Real) = RK4ErrorValue(_to_interval(x), 0.0)
Base.convert(::Type{RK4ErrorValue}, x::Real) = _error_value(x)
Base.promote_rule(::Type{RK4ErrorValue}, ::Type{Float64}) = RK4ErrorValue

function _error_value(range, error, rounded_operands_range)
    total = error + _roundoff_radius(rounded_operands_range)
    return RK4ErrorValue(range, _valid_interval_enclosure(total) ? IA.sup(total) : Inf)
end

function Base.:+(a::RK4ErrorValue, b::RK4ErrorValue)
    range = a.range + b.range
    isfinite(a.error) && isfinite(b.error) || return RK4ErrorValue(range, Inf)
    return _error_value(
        range,
        _to_interval(a.error) + _to_interval(b.error),
        _expanded(a) + _expanded(b),
    )
end
function Base.:-(a::RK4ErrorValue, b::RK4ErrorValue)
    range = a.range - b.range
    isfinite(a.error) && isfinite(b.error) || return RK4ErrorValue(range, Inf)
    return _error_value(
        range,
        _to_interval(a.error) + _to_interval(b.error),
        _expanded(a) - _expanded(b),
    )
end
Base.:-(a::RK4ErrorValue) = RK4ErrorValue(-a.range, a.error)
function Base.:*(a::RK4ErrorValue, b::RK4ErrorValue)
    range = a.range * b.range
    isfinite(a.error) && isfinite(b.error) || return RK4ErrorValue(range, Inf)
    ea, eb = _to_interval(a.error), _to_interval(b.error)
    error = abs(a.range) * eb + abs(b.range) * ea + ea * eb
    return _error_value(range, error, _expanded(a) * _expanded(b))
end
function Base.:/(a::RK4ErrorValue, b::RK4ErrorValue)
    range = a.range / b.range
    isfinite(a.error) && isfinite(b.error) || return RK4ErrorValue(range, Inf)
    divisor = _expanded(b)
    error = (_to_interval(a.error) + abs(range) * b.error) / abs(divisor)
    return _error_value(range, error, _expanded(a) / divisor)
end
Base.abs(a::RK4ErrorValue) = RK4ErrorValue(abs(a.range), a.error)

for f in (:+, :-, :*, :/)
    @eval begin
        Base.$f(a::RK4ErrorValue, b::Real) = $f(a, _error_value(b))
        Base.$f(a::Real, b::RK4ErrorValue) = $f(_error_value(a), b)
    end
end

# Elementary functions use a derivative bound on the union of real and
# rounded arguments. If the derivative cannot be enclosed, the direct
# interval RK4 image remains available.
_expanded(a::RK4ErrorValue) = a.range + IA.interval(Float64, -a.error, a.error)

function _elementary_derivative(f, x)
    f === sin && return cos(x)
    f === cos && return -sin(x)
    f === tan && return 1 / cos(x)^2
    f === exp && return exp(x)
    f === expm1 && return exp(x)
    f === log && return 1 / x
    f === log1p && return 1 / (1 + x)
    f === sinh && return cosh(x)
    f === cosh && return sinh(x)
    f === tanh && return 1 - tanh(x)^2
    f === asin && return 1 / sqrt(1 - x^2)
    f === acos && return -1 / sqrt(1 - x^2)
    f === atan && return 1 / (1 + x^2)
    return IA.entireinterval(Float64)
end

function (op::RoundedElementary)(a::RK4ErrorValue)
    range = op(a.range)
    isfinite(a.error) || return RK4ErrorValue(range, Inf)
    x = _expanded(a)
    error = abs(_elementary_derivative(op.f, x)) * a.error
    return _error_value(range, error, op(x))
end

function _power_derivative(op, x, n)
    iszero(n) && return _to_interval(0.0)
    # Subtract the exponent in interval arithmetic, unless it is a small
    # integer whose predecessor can be represented exactly.
    power = if isinteger(n) && abs(n) < 2.0^52
        op(x, Int(n) - 1)
    else
        IA.pow(x, _to_interval(n) - _to_interval(1.0))
    end
    return _to_interval(n) * power
end

function (op::RoundedElementary{typeof(^)})(a::RK4ErrorValue, n::Real)
    range = op(a.range, n)
    isfinite(a.error) || return RK4ErrorValue(range, Inf)
    x = _expanded(a)
    error = abs(_power_derivative(op, x, n)) * a.error
    return _error_value(range, error, op(x, n))
end

function Base.sqrt(a::RK4ErrorValue)
    range = sqrt(a.range)
    isfinite(a.error) || return RK4ErrorValue(range, Inf)
    x = _expanded(a)
    rounded_range = sqrt(x)
    return _error_value(range, a.error / (2 * rounded_range), rounded_range)
end

struct RK4Derivative <: Real
    value::IA.Interval{Float64}
    gradient::Vector{IA.Interval{Float64}}
end

function Base.:+(a::RK4Derivative, b::RK4Derivative)
    return RK4Derivative(a.value + b.value, a.gradient .+ b.gradient)
end
function Base.:-(a::RK4Derivative, b::RK4Derivative)
    return RK4Derivative(a.value - b.value, a.gradient .- b.gradient)
end
Base.:-(a::RK4Derivative) = RK4Derivative(-a.value, .-a.gradient)
function Base.:*(a::RK4Derivative, b::RK4Derivative)
    return RK4Derivative(a.value * b.value, a.gradient .* b.value .+ a.value .* b.gradient)
end
function Base.:/(a::RK4Derivative, b::RK4Derivative)
    ratio = a.value / b.value
    return RK4Derivative(ratio, (a.gradient .- ratio .* b.gradient) ./ b.value)
end

for f in (:+, :-, :*, :/)
    @eval begin
        function Base.$f(a::RK4Derivative, b::Real)
            return $f(a, RK4Derivative(_to_interval(b), zero(a.gradient)))
        end
        function Base.$f(a::Real, b::RK4Derivative)
            return $f(RK4Derivative(_to_interval(a), zero(b.gradient)), b)
        end
    end
end

Base.abs(a::RK4Derivative) =
    RK4Derivative(abs(a.value), IA.interval(Float64, -1, 1) .* a.gradient)
Base.sqrt(a::RK4Derivative) =
    RK4Derivative(sqrt(a.value), a.gradient ./ (2 * sqrt(a.value)))

function (op::RoundedElementary)(a::RK4Derivative)
    return RK4Derivative(op(a.value), _elementary_derivative(op.f, a.value) .* a.gradient)
end

function (op::RoundedElementary{typeof(^)})(a::RK4Derivative, n::Real)
    return RK4Derivative(op(a.value, n), _power_derivative(op, a.value, n) .* a.gradient)
end

const _rk4_elementary =
    (sin, cos, tan, exp, expm1, log, log1p, sinh, cosh, tanh, asin, acos, atan)

for f in _rk4_elementary
    correct = getfield(CoreMath, Symbol(:cr_, nameof(f)))
    @eval function (op::RoundedElementary{typeof($f)})(x::Float64)
        return $correct(x)
    end
end

(::RoundedElementary{typeof(^)})(x::Float64, y::Float64) = CoreMath.cr_pow(x, y)
(op::RoundedElementary{typeof(^)})(x::Float64, n::Int) = op(x, Float64(n))

(op::RoundedElementary)(x::IA.Interval{Float64}) = op.f(x)

(::RoundedElementary{typeof(^)})(x::IA.Interval{Float64}, n::Integer) = IA.pown(x, n)
(op::RoundedElementary{typeof(^)})(x::IA.Interval{Float64}, n::Float64) =
    isinteger(n) ? op(x, Integer(n)) : IA.pow(x, _to_interval(n))

# With `cse = true`, Symbolics hoists every literal into a variable that is bound once,
# e.g. `var"##cse#1" = 2` and then `x^var"##cse#1"`. Record the literals bound exactly
# once, so that an exponent hoisted this way is still seen as a constant.
function _literal_bindings(ex)
    counts = Dict{Symbol, Int}()
    literals = Dict{Symbol, Float64}()
    function visit(e)
        e isa Expr || return nothing
        if e.head == :(=) && e.args[1] isa Symbol
            name = e.args[1]
            counts[name] = get(counts, name, 0) + 1
            e.args[2] isa Real && (literals[name] = Float64(e.args[2]))
        end
        foreach(visit, e.args)
        return nothing
    end
    visit(ex)
    return Dict(name => value for (name, value) in literals if counts[name] == 1)
end

# Symbolics emits resolved function objects in arithmetic calls. Keep the
# generated bindings and array construction, but reject unknown operations.
_certified_rk4_expression(ex) = _certified_rk4_expression(ex, _literal_bindings(ex))

function _certified_rk4_expression(ex, literals)
    if ex isa Real
        value = Float64(ex)
        isfinite(value) || throw(ArgumentError("nonfinite coefficient in certified RK4"))
        return value
    end
    ex isa Expr || return ex
    ex.head == :ref && return ex # generated argument indexing
    ex.head in (:if, :ifelse, :&&, :||, :comparison) && throw(
        ArgumentError("certified RK4 dynamics must be a branch-free symbolic expression"),
    )
    convert_all(args) = map(a -> _certified_rk4_expression(a, literals), args)
    if ex.head == :call
        op = ex.args[1]
        op === getindex && return ex
        if op === Symbolics.SymbolicUtils.Code.create_array
            return Expr(:call, ex.args[1:5]..., convert_all(ex.args[6:end])...)
        end
        args = convert_all(ex.args[2:end])
        if op === (^)
            # The exponent is a numeric literal or a variable bound once to one.
            length(args) == 2 ||
                throw(ArgumentError("certified RK4 requires a binary power"))
            exponent = args[2] isa Symbol ? get(literals, args[2], nothing) : args[2]
            exponent isa Float64 ||
                throw(ArgumentError("certified RK4 requires a constant exponent"))
            return Expr(:call, _rk4_operation, RoundedElementary(op), args[1], exponent)
        elseif op in _rk4_elementary
            length(args) == 1 || throw(ArgumentError("certified RK4 requires unary $op"))
            return Expr(:call, _rk4_operation, RoundedElementary(op), args...)
        end
        op in (+, -, *, /, abs, sqrt, getindex) ||
            throw(ArgumentError("unsupported operation in certified RK4 dynamics: $op"))
        # Lower n-ary arithmetic to the same binary operation graph for all
        # scalar types. Julia's generic n-ary methods need not group identically.
        if op in (+, *) && length(args) > 2
            return foldl((a, b) -> Expr(:call, _rk4_operation, op, a, b), args)
        end
        return Expr(:call, _rk4_operation, op, args...)
    end
    return Expr(ex.head, convert_all(ex.args)...)
end

function _rk4_step_jacobian(f, stages, u, h)
    n = length(first(stages))
    identity = [_to_interval(i == j ? 1.0 : 0.0) for i in 1:n, j in 1:n]
    input = [RK4Derivative(v, fill(_to_interval(0.0), n)) for v in u]
    function jacobian_at(y)
        state = [RK4Derivative(y[i], identity[i, :]) for i in 1:n]
        result = f(state, input)
        return [
            result[i] isa RK4Derivative ? result[i].gradient[j] : _to_interval(0.0) for
            i in 1:n, j in 1:n
        ]
    end
    J1, J2, J3, J4 = map(jacobian_at, stages)
    all(J -> all(_valid_interval_enclosure, J), (J1, J2, J3, J4)) || return nothing
    product(A, B) =
        [sum(A[i, k] * B[k, j] for k in axes(A, 2)) for i in axes(A, 1), j in axes(B, 2)]
    d2 = product(J2, identity .+ (h / 2) .* J1)
    d3 = product(J3, identity .+ (h / 2) .* d2)
    d4 = product(J4, identity .+ h .* d3)
    return identity .+ (J1 .+ 2 .* d2 .+ 2 .* d3 .+ d4) .* (h / 6)
end

function _centered_rk4_image(approx, source, stages, ubox, h, f)
    derivative = _rk4_step_jacobian(approx.f, stages, ubox, h)
    derivative === nothing && return nothing
    all(_valid_interval_enclosure, derivative) || return nothing
    center = _point_box(IA.mid.(source))
    result = _try_rk4_stages(f, center, ubox, h)
    result === nothing && return nothing
    _, _, center_image = result

    function error_f(x, u)
        return _error_value.(approx.f(x, u))
    end
    # Use the same stored Float64 coefficients as the concrete map.
    # Arithmetic on RK4ErrorValue encloses each roundoff.
    values = _try_rk4_stages(error_f, _error_value.(source), _error_value.(ubox), h)
    values === nothing && return nothing
    _, _, values = values
    all(v -> isfinite(v.error), values) || return nothing
    displacement = source .- center
    return map(eachindex(center_image)) do i
        error = values[i].error
        return center_image[i] +
               sum(derivative[i, j] * displacement[j] for j in eachindex(source)) +
               IA.interval(Float64, -error, error)
    end
end

function _compile_certified_rk4(f, x, u)
    convert_expression(ex, state) = _certified_rk4_expression(Symbolics.toexpr(ex, state))
    return first(
        Symbolics.build_function(
            f,
            x,
            u;
            expression = Val(false),
            nanmath = false,
            cse = true,
            conv = convert_expression,
        ),
    )
end

# The endpoint must retain its outward-rounded bounds. Converting to a
# center-radius Hyperrectangle can move a face inward.
struct CertifiedRK4Image{V} <: LazySets.LazySet{Float64}
    low::V
    high::V
end

LazySets.low(S::CertifiedRK4Image) = S.low
LazySets.high(S::CertifiedRK4Image) = S.high
LazySets.dim(S::CertifiedRK4Image) = length(S.low)

function _certified_grid(mapping)
    mapping isa Union{MP.ExplicitGridMapping, MP.ImplicitGridMapping} || throw(
        ArgumentError(
            "certified RK4 requires a nonperiodic explicit or implicit grid mapping",
        ),
    )
    return _certified_grid_object(MP.get_grid(mapping))
end

function _certified_grid_object(grid)
    grid isa MP.GridFree &&
    eltype(MP.get_h(grid)) === Float64 &&
    all(isfinite, MP.get_origin(grid)) &&
    all(h -> isfinite(h) && h > 0, MP.get_h(grid)) || throw(
        ArgumentError(
            "certified RK4 requires a Float64 GridFree with positive finite spacing",
        ),
    )
    return grid
end

function _quantized_positions(grid, image::CertifiedRK4Image)
    # Subtraction, division by positive spacing, and ties-to-even rounding
    # are all monotone. Evaluate the very same Float64 quantizer at both ends.
    local lo, hi
    try
        lo = MP.get_pos_by_coord(grid, image.low)
        hi = MP.get_pos_by_coord(grid, image.high)
    catch err
        err isa InexactError || rethrow()
        return nothing
    end
    return map((a, b) -> a:b, lo, hi)
end

function MP.get_states_from_set_strict(
    mapping::MP.GridMapping,
    image::CertifiedRK4Image,
    mode::MP.INCL_MODE,
)
    mode == MP.OUTER || throw(ArgumentError("certified RK4 images support OUTER lookup"))
    positions = _quantized_positions(_certified_grid(mapping), image)
    targets = Int[]
    positions === nothing && return targets, false
    for pos in Iterators.product(positions...)
        MP.is_valid_pos(mapping, pos) || return targets, false
        push!(targets, MP.get_state_by_pos(mapping, pos))
    end
    return targets, true
end

# Invert the three rounded operations of the Float64 quantizer conservatively.
# Padding a range by one representable value encloses all exact operands that
# could have rounded into the range, including at zero and binade boundaries.
function _rounding_preimage(a)
    return IA.interval(Float64, prevfloat(IA.inf(a)), nextfloat(IA.sup(a)))
end

# Order Float64 bit patterns as unsigned integers so bisection also works
# across zero and at the ends of the finite range.
const _float_sign = UInt64(1) << 63
function _float_key(x::Float64)
    bits = reinterpret(UInt64, x)
    return signbit(x) ? ~bits : bits ⊻ _float_sign
end

function _key_float(key::UInt64)
    bits = iszero(key & _float_sign) ? ~key : key ⊻ _float_sign
    return reinterpret(Float64, bits)
end

function _first_float_satisfying(test, left, right)
    test(left) && return left
    test(right) || return nothing
    lo, hi = _float_key(left), _float_key(right)
    while hi - lo > 1
        mid = lo + (hi - lo) ÷ 2
        if test(_key_float(mid))
            hi = mid
        else
            lo = mid
        end
    end
    return _key_float(hi)
end

function _quantizer_coordinate(pos, origin, spacing)
    half = IA.interval(Float64, -0.5, 0.5)
    quotient_input = _rounding_preimage(_to_interval(pos) + half)
    subtraction_input = _rounding_preimage(quotient_input * _to_interval(spacing))
    outer = subtraction_input + _to_interval(origin)
    _valid_interval_enclosure(outer) || return outer

    # The Float64 result of round has the same value as round(Int, ...)
    # whenever the latter is defined. It can also compare out-of-range results.
    quantize(x) = round((x - origin) / spacing)
    lo, hi = IA.inf(outer), IA.sup(outer)
    first = _first_float_satisfying(x -> quantize(x) >= pos, lo, hi)
    first === nothing && return IA.emptyinterval(Float64)
    quantize(first) == pos || return IA.emptyinterval(Float64)
    quantize(hi) <= pos && return IA.interval(Float64, first, hi)
    after = _first_float_satisfying(x -> quantize(x) > pos, first, hi)
    return IA.interval(Float64, first, prevfloat(after))
end

function _quantizer_cell(grid, pos)
    origin, spacing = MP.get_origin(grid), MP.get_h(grid)
    return [_quantizer_coordinate(pos[i], origin[i], spacing[i]) for i in eachindex(pos)]
end

# Singleton slope enclosures mean every concrete execution sees exactly these
# Float64 slopes. The final update in this coordinate is then monotone in x.
function _constant_slope_image(source, slopes, h, direct)
    return map(eachindex(source)) do i
        ks = map(k -> k[i], slopes)
        all(IA.isthin, ks) || return direct[i]
        k1, k2, k3, k4 = map(IA.inf, ks)
        increment = (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6)
        return IA.interval(
            Float64,
            IA.inf(source[i]) + increment,
            IA.sup(source[i]) + increment,
        )
    end
end

function CertifiedRK4GrowthBound(
    f::AbstractVector,
    x::AbstractVector,
    u::AbstractVector,
    X,
    U,
    tstep::Real;
    nsystem::Int = ST.DEFAULT_NUM_SUBSTEPS,
)
    pkgversion(IA) >= v"1" ||
        throw(ArgumentError("certified RK4 requires IntervalArithmetic 1.x"))
    length(f) == length(x) ||
        throw(DimensionMismatch("one expression is required per state"))
    τ = Float64(tstep)
    nsystem > 0 && isfinite(τ) && τ >= 0 || throw(
        ArgumentError(
            "certified RK4 requires a finite nonnegative step and positive substep count",
        ),
    )
    fsym, xs, us = Symbolics.Num.(f), Symbolics.Num.(x), Symbolics.Num.(u)
    compiled = _compile_certified_rk4(fsym, xs, us)
    # Identifies the dynamics in the provenance of the transitions (`get_provenance`): the printed
    # expressions are the same in every session, unlike the compiled function.
    dynamics_fingerprint =
        fingerprint_string(fingerprint((string.(fsym), string.(xs), string.(us))))
    n = length(x)
    m = length(us)
    dynamics =
        (x, u) -> SVector{n, Float64}(compiled(Vector{Float64}(x), Vector{Float64}(u)))
    h = τ / nsystem
    function map_float64(x, u)
        x, u = SVector{n, Float64}(x), SVector{m, Float64}(u)
        for _ in 1:nsystem
            _, _, x = _rk4_stages(dynamics, x, u, h)
        end
        return x
    end
    system = MS.ConstrainedBlackBoxControlDiscreteSystem(map_float64, n, m, X, U)
    return CertifiedRK4GrowthBound(system, compiled, τ, nsystem, dynamics_fingerprint)
end

function CertifiedRK4GrowthBound(
    system::MS.ConstrainedBlackBoxControlContinuousSystem,
    tstep::Real;
    kwargs...,
)
    x = [Symbolics.variable(:x, i) for i in 1:MS.statedim(system)]
    u = [Symbolics.variable(:u, i) for i in 1:MS.inputdim(system)]
    f = collect(MS.mapping(system)(x, u))
    return CertifiedRK4GrowthBound(
        f,
        x,
        u,
        MS.stateset(system),
        MS.inputset(system),
        tstep;
        kwargs...,
    )
end

function _certified_rk4_image(approx, source, ubox)
    IA.default_rounding() isa IA.IntervalRounding{:correct} ||
        throw(ArgumentError("certified RK4 requires correct interval rounding"))
    rounding(Float64) == RoundNearest && !get_zero_subnormals() || throw(
        ArgumentError("certified RK4 requires round-to-nearest with subnormals enabled"),
    )
    all(_valid_interval_enclosure, source) && all(_valid_interval_enclosure, ubox) ||
        return nothing
    function f(x, u)
        return _checked_primitive.(_to_interval.(approx.f(x, u)))
    end
    h = approx.tstep / approx.nsystem
    endpoint = source
    for _ in 1:approx.nsystem
        result = _try_rk4_stages(f, endpoint, ubox, h)
        result === nothing && return nothing
        stages, slopes, direct = result
        centered = _centered_rk4_image(approx, endpoint, stages, ubox, h, f)
        direct = _constant_slope_image(endpoint, slopes, h, direct)
        endpoint = if centered === nothing || !all(_valid_interval_enclosure, centered)
            direct
        else
            [_box_intersection(a, b) for (a, b) in zip(direct, centered)]
        end
        all(_valid_interval_enclosure, endpoint) || return nothing
    end
    return CertifiedRK4Image(IA.inf.(endpoint), IA.sup.(endpoint))
end

_box_intersection(a, b) =
    IA.interval(Float64, max(IA.inf(a), IA.inf(b)), min(IA.sup(a), IA.sup(b)))

# Image of the box `elem` itself. A Float64 state that the grid quantizer assigns to a cell can lie a
# few ulps outside the nominal box of that cell, so this is not the image to use for a grid
# abstraction: see `overapproximation_map`.
function ST.reach_set(approx::CertifiedRK4GrowthBound, elem, u, cache)
    center, radius = LazySets.center(elem), LazySets.radius_hyperrectangle(elem)
    source = [
        IA.interval(Float64, center[i], radius[i]; format = :midpoint) for
        i in 1:length(center)
    ]
    return _certified_rk4_image(approx, source, _point_box(Float64.(u)))
end

ST.get_over_approximation_map(approx::CertifiedRK4GrowthBound) =
    (elem, u) -> ST.reach_set(approx, elem, u, nothing)

# The action is disabled: no transition at all for this pair. Dionysos recognizes it by the lookup of
# the target cells (`allin = false`), which every version of the kernel honors. (Recent versions
# also accept `nothing` from a reach set, older ones do not.)
struct DisabledImage <: LazySets.LazySet{Float64}
    dim::Int
end

LazySets.dim(S::DisabledImage) = S.dim

function MP.get_states_from_set_strict(::MP.GridMapping, ::DisabledImage, ::MP.INCL_MODE)
    return Int[], false
end

"""
    overapproximation_map(approx::CertifiedRK4GrowthBound, state_grid) -> Function

The function `(cell, u) -> image` to give to the `UniformGridAbstraction` optimizer as its
`overapproximation_map` (with `approx_mode = USER_DEFINED`). `state_grid` is the nonperiodic Float64
`GridFree` of the abstraction. The source set is not the nominal box of `cell` but the exact set of
Float64 states that the grid quantizer assigns to it, and the image is looked up with the same
quantizer, ties to even included. An action that cannot be certified inside the mapped domain gets
no transition.
"""
function overapproximation_map(approx::CertifiedRK4GrowthBound, state_grid)
    _certified_grid_object(state_grid)
    N = MP.get_dim(state_grid)
    return (elem, u) -> begin
        pos = MP.get_pos_by_coord(state_grid, LazySets.center(elem))
        source = _quantizer_cell(state_grid, pos)
        image = _certified_rk4_image(approx, source, _point_box(Float64.(u)))
        return image === nothing ? DisabledImage(N) : image
    end
end

end
