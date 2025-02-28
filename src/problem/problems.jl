"""
The structure 

    ProblemType

is the abstract type which defines a problem.
"""
abstract type ProblemType end

"""
The structure 

    EmptyProblem{S, X}

encodes a problem consisting in constructing the abstraction in X 
- `S` is the system,
- `X` is the region, 
"""
struct EmptyProblem{S, X} <: ProblemType
    system::S
    region::X
end

"""
The structure 

    OptimalControlProblem{S, XI, XT, XC, TC, T}

encodes an optimal control problem where 
- `S` is the system,
- `XI` is the initial set, 
- `XT` is the target set,
- `XC` is the state cost,
- `TC` is transistion cost and
- `T` is the number of allowed time steps
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
- `S` is the system,
- `XI` is the initial set, 
- `XS` is the safe set and
- `T` is the number of allowed time steps
"""
struct SafetyProblem{S, XI, XS, T <: Real} <: ProblemType
    system::S
    initial_set::XI
    safe_set::XS
    time::T
end

struct Infinity <: Real end
Base.isfinite(::Infinity) = false

@recipe function f(problem::EmptyProblem; domain_color = :gray, region_color = :lightgray)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Region"
        color := region_color
        problem.region
    end
end

@recipe function f(
    problem::OptimalControlProblem;
    domain_color = :gray,
    initial_set_color = :green,
    target_set_color = :red,
)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
    @series begin
        label := "Target set"
        color := target_set_color
        problem.target_set
    end
end

@recipe function f(
    problem::SafetyProblem;
    domain_color = :gray,
    initial_set_color = :green,
    safe_set_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        problem.system.X
    end
    @series begin
        label := "Safe set"
        color := safe_set_color
        problem.safe_set
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
end

export OptimalControlProblem
export SafetyProblem
export Infinity
