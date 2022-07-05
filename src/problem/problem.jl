module Problem


    struct Infinity <: Real end


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

    include("optim/bemporad_morari.jl")
    include("optim/branch_and_bound.jl")
    include("optim/q_learning.jl")

    export OptimalControlProblem
    export SafetyProblem
    export Infinity
end
