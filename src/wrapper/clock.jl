# ----------------------------------------------------------------------------------------
# Clocks.
#
# A clock is not a new kind of declaration: it is a state whose dynamics is `∂(t) == 1` where
# it runs and `∂(t) == 0` where it is frozen. That is exactly what
# `Symbolic.ClockAbstraction` accepts (`A = I` or `A = 0`), so the DSL needs no new syntax —
# only the recognition rule below.
#
# A model with a clock has augmented states `(x, t, mode)`: the clock is a state of the hybrid
# system but not a coordinate of the physical `x`, so it is kept out of `state_indices`.
# ----------------------------------------------------------------------------------------

# The constant rate of a clock equation, or `nothing` when the right-hand side depends on a
# variable (and is therefore ordinary dynamics).
#
# The right-hand side is not a bare number by the time it reaches us: JuMP wraps a literal, so
# `∂(t) == 1` arrives as the expression `+(1.0)`.
_clock_rate(rate::Number) = Float64(rate)
_clock_rate(::Any) = nothing

function _clock_rate(f::MOI.ScalarNonlinearFunction)
    f.head === :+ || return nothing
    total = 0.0
    for arg in f.args
        value = _clock_rate(arg)
        value === nothing && return nothing
        total += value
    end
    return total
end

function _clock_rate(f::MOI.ScalarAffineFunction)
    isempty(f.terms) || return nothing
    return f.constant
end

# Every rate declared for variable `i`, across the model level and all modes. `nothing` if any
# of them is not constant.
function _declared_rates(ir::ModelIR, i::Int)
    rates = Float64[]
    d = ir.dynamics[i]
    if d !== nothing
        rate = _clock_rate(d)
        rate === nothing && return nothing
        push!(rates, rate)
    end
    for m in values(ir.modes)
        haskey(m.dynamics, i) || continue
        rate = _clock_rate(m.dynamics[i])
        rate === nothing && return nothing
        push!(rates, rate)
    end
    return rates
end

"""
    detect_clock!(ir::ModelIR) -> Union{Nothing, Int}

Find the clock (rule I2) and mark it `CLOCK`, returning its variable index.

A variable is a clock when every equation declared for it is the constant `1` (running) or `0`
(frozen), at least one of them is `1` — a clock has to run somewhere — and it drives no other
state. That last condition is what separates a clock from an ordinary state that merely happens
to be held constant.
"""
function detect_clock!(ir::ModelIR)
    is_hybrid(ir) || return nothing

    # Variables appearing on some *other* variable's right-hand side: those are inputs to the
    # physical dynamics, not clocks.
    driving = Set{Int}()
    for d in ir.dynamics
        d === nothing || _collect_variables!(driving, d)
    end
    for m in values(ir.modes), d in values(m.dynamics)
        _collect_variables!(driving, d)
    end

    candidates = Int[]
    for i in eachindex(ir.variables)
        i in driving && continue
        rates = _declared_rates(ir, i)
        rates === nothing && continue
        isempty(rates) && continue
        all(r -> r == 0.0 || r == 1.0, rates) || continue
        any(==(1.0), rates) || continue
        push!(candidates, i)
    end

    isempty(candidates) && return nothing
    length(candidates) == 1 || error(
        "Several variables look like clocks ($(join((describe(ir.variables[i], i) for i in candidates), ", "))); " *
        "a hybrid model carries at most one time axis.",
    )

    index = candidates[]
    ir.variables[index].role = CLOCK
    return index
end

"""
    clock_index(ir) -> Union{Nothing, Int}

The clock variable, or `nothing` for a time-free model.
"""
function clock_index(ir::ModelIR)
    index = findfirst(v -> v.role === CLOCK, ir.variables)
    return index
end

"""
    clock_system(ir, mode, index) -> MathematicalSystems.ConstrainedLinearContinuousSystem

The clock subsystem of one mode: `ẋ = 1` where the clock runs, `ẋ = 0` where the mode freezes
it, over the time domain given by the clock variable's bounds.
"""
function clock_system(ir::ModelIR, m::ModeIR, index::Int)
    expr = get(m.dynamics, index, ir.dynamics[index])
    rate = expr === nothing ? nothing : _clock_rate(expr)
    rate === nothing && error(
        "Mode $(m.id) does not say whether the clock $(describe(ir.variables[index], index)) " *
        "runs: add `∂(t) == 1` or `∂(t) == 0` to it.",
    )
    lo, hi = _mode_bound(ir, m, index)
    domain = UT.box(SVector(lo), SVector(hi))
    return MS.ConstrainedLinearContinuousSystem(SMatrix{1, 1}(rate), domain)
end
