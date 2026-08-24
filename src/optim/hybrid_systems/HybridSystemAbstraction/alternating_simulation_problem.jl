
mutable struct OptimizerAlternatingSimulationProblem{T} <: OP.AbstractDionysosOptimizer
    # Input
    alternating_simulation_problem::Union{Nothing, PR.AlternatingSimulationProblem}
    concrete_system::Union{Nothing, HybridSystem}
    optimizer_list::Union{Nothing, AbstractVector}
    optimizer_kwargs_dict::Union{Nothing, AbstractVector}
    # Per-mode `nothing` (build) or the index of an earlier mode whose abstraction
    # this mode reuses — see `build_mode_symbolic_models`.
    shared_abstraction::Union{Nothing, AbstractVector}
    # Abstract the distinct modes on separate threads. Opt-in: nesting this inside a mode's own
    # threaded build backend oversubscribes rather than accelerates.
    parallel_modes::Bool

    max_iterations::Union{Nothing, Int}
    print_level::Int

    # Output 
    abstract_system::Union{Nothing, SY.HybridSymbolicModel}
    abstraction_construction_time_sec::T

    function OptimizerAlternatingSimulationProblem{T}() where {T}
        optimizer = new{T}(
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            1000,
            1,
            nothing,
            0.0,
        )
        return optimizer
    end
end

OptimizerAlternatingSimulationProblem() = OptimizerAlternatingSimulationProblem{Float64}()

MOI.is_empty(optimizer::OptimizerAlternatingSimulationProblem) =
    optimizer.alternating_simulation_problem === nothing

function MOI.get(model::OptimizerAlternatingSimulationProblem, ::MOI.SolveTimeSec)
    return model.abstraction_construction_time_sec
end

function reset!(model::OptimizerAlternatingSimulationProblem)
    return model
end

function _validate_model(
    model::OptimizerAlternatingSimulationProblem,
    required_fields::Vector{Symbol},
)
    for field in required_fields
        if isnothing(getfield(model, field))
            error(
                "Please set the `$(field)`. Missing required field in OptimizerAlternatingSimulationProblem.",
            )
        end
    end
end

# Losing part of a switch is legitimate — the guard cells at the edge of the target mode's
# domain have nowhere to land — but it is the kind of thing that explains an unexpectedly small
# winning set an hour later, so it is summarised rather than left in the log.
function _print_build_report(report::SY.HybridBuildReport)
    isempty(report) && return
    for (id, (dropped, total)) in sort!(collect(report.dropped_resets))
        println("  transition $id: $dropped of $total guard cells reset out of domain")
    end
    for (id, offset) in sort!(collect(report.inexact_resets))
        println(
            "  transition $id: reset is not lattice-exact (offset ",
            round(offset; sigdigits = 3),
            ") — the abstraction may be unsound",
        )
    end
    return
end

function MOI.optimize!(optimizer::OptimizerAlternatingSimulationProblem)
    t_ref = time()
    # Ensure necessary parameters are set
    _validate_model(optimizer, [:alternating_simulation_problem])
    @assert optimizer.alternating_simulation_problem.system !== nothing "System must be set to construct the abstraction."

    # Build symbolic model
    optimizer.print_level >= 1 &&
        println("Construct the Hybrid System Abstraction: started")
    optimizer.abstract_system = build_timed_hybrid_symbolic_model(
        optimizer.alternating_simulation_problem.system,
        optimizer.optimizer_list,
        optimizer.optimizer_kwargs_dict;
        shared_abstraction = optimizer.shared_abstraction,
        parallel_modes = optimizer.parallel_modes,
    )

    optimizer.print_level >= 1 && println(
        "Construct the Hybrid System Abstraction: terminated with success: ",
        "$(HybridSystems.ntransitions(SY.get_automaton(optimizer.abstract_system))) transitions created",
    )
    optimizer.print_level >= 1 &&
        _print_build_report(SY.get_build_report(optimizer.abstract_system))

    optimizer.abstraction_construction_time_sec = time() - t_ref
    return
end
