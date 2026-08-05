# ----------------------------------------------------------------------------------------
# Which solver runs the lowered problem.
#
# The choice follows from the shape of the model, and is overridable. Both hooks are ordinary
# multiple dispatch, so a new solver family joins by adding one method to each — the same
# pattern as `Optim.control_solver_for`.
# ----------------------------------------------------------------------------------------

const UniformGrid = OP.Abstraction.UniformGridAbstraction.Optimizer
const HybridAbstraction = OP.Abstraction.HybridSystemAbstraction.Optimizer

"""
    select_solver(system, problem) -> optimizer type

The solver family to use for a lowered `(system, problem)` pair when the user set no
`"solver"` attribute. Add a method to support a new family.
"""
select_solver(::Any, ::Any) = UniformGrid
select_solver(::HybridSystems.HybridSystem, ::Any) = HybridAbstraction

"""
    supports_problem(solver_type, problem_type) -> Bool

Whether `solver_type` can solve `problem_type`. Declared per family so an unsupported
combination is reported by the front-end, naming both, instead of failing deep inside a
sub-solver.
"""
supports_problem(::Type, ::Type{<:PR.ProblemType}) = false

function supports_problem(
    ::Type{<:UniformGrid},
    ::Type{
        <:Union{
            PR.OptimalControlProblem,
            PR.SafetyProblem,
            PR.ReachAndStayProblem,
            PR.CoSafeLTLProblem,
            PR.AlternatingSimulationProblem,
        },
    },
)
    return true
end

function supports_problem(
    ::Type{<:HybridAbstraction},
    ::Type{<:Union{PR.OptimalControlProblem, PR.SafetyProblem}},
)
    return true
end

# The solver type the model should run on, given the problem it lowered to.
function _solver_type(model, problem)
    factory = model.solver_factory
    factory === nothing && return select_solver(problem.system, problem)
    return factory
end

function _check_supported(solver_type, problem)
    supports_problem(solver_type, typeof(problem)) && return
    return error(
        "$(solver_type) does not support $(typeof(problem).name.name). Choose another " *
        "solver with `set_attribute(model, \"solver\", …)`, or state a different " *
        "specification.",
    )
end
