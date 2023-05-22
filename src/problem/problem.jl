module Problem

struct Infinity <: Real end

Base.isfinite(::Infinity) = false

struct OptimalControlProblem{S, XI, XT, XC, TC, T<:Real}
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
end

struct SafetyProblem{S, XI, XS, T<:Real}
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

export OptimalControlProblem
export SafetyProblem
export Infinity

end
