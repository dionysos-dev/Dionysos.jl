export ZeroCost, ConstantCost, QuadraticControlCost
export OptimalControlProblem

struct ZeroCost end

struct ConstantCost{T}
    cost::T
end

struct QuadraticControlCost{T, MT <: AbstractMatrix{T}}
    Q::MT
end

struct OptimalControlProblem{S, Q0, X0, XC, TC}
    system::S
    q_0::Q0
    x_0::X0
    state_cost::XC
    transition_cost::TC
    number_of_time_steps::Int
end

export optimal_control
function optimal_control end
