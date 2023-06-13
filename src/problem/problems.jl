abstract type ProblemType end

struct OptimalControlProblem{S, XI, XT, XC, TC, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
end

struct SafetyProblem{S, XI, XS, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

struct Infinity <: Real end
Base.isfinite(::Infinity) = false

export OptimalControlProblem
export SafetyProblem
export Infinity


