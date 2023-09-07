abstract type ProblemType end


"""
The structure 

    OptimalControlProblem{S, XI, XT, XC, TC, T}

encodes an optimal control problem where 
`S` is the system,
`XI` is the initial set, 
`XT` is the target set,
`XC` is the state cost,
`TC` is transistion cost and
`T` is the number of allowed time steps
"""
struct OptimalControlProblem{S, XI, XT, XC, TC, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    target_set::XT
    state_cost::XC
    transition_cost::TC
    time::T
end

"""
The structure 

    SafetyProblem{S, XI, XS, T}

encodes a safety problem where
`S` is the system,
`XI` is the initial set, 
`XS` is the safe set and
`T` is the number of allowed time steps
"""
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
