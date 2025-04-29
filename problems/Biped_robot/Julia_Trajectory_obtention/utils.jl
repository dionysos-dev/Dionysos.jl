using MathematicalSystems, StaticArrays, Plots, LinearAlgebra
using JuMP
using JLD2

# include Dionysos
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic
const OP = DI.Optim
const AB = OP.Abstraction

function get_abstract_closed_loop_trajectory(abstract_system, abstract_controller, source, nstep; stopping = (s)->false)
    state_traj, input_traj = [source], []

    for _ in 1:nstep
        stopping(source) && break
        abstract_inputs = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(abstract_inputs)
            @warn("Uncontrollable state: $source")
            return nothing
        end
        input = first(abstract_inputs)[1]
        targets = []
        UT.fix_and_eliminate_tail!(targets, abstract_system.autom.transitions, (source, input))
        source = first(targets)
        push!(state_traj, source)
        push!(input_traj, input)
    end
    return ST.Control_trajectory(ST.Trajectory(state_traj), ST.Trajectory(input_traj))
end

function get_concrete_trajectory(abstract_system, abstract_trajectory::ST.Control_trajectory)
    concrete_state_traj, concrete_input_traj = [], []
    for k in 1:ST.length(abstract_trajectory)
        abstarct_state = ST.get_state(abstract_trajectory, k)
        concrete_state = SY.get_concrete_state(abstract_system, abstarct_state)
        push!(concrete_state_traj, concrete_state)

        if k<ST.length(abstract_trajectory)
            abstract_input = ST.get_input(abstract_trajectory, k)
            concrete_input = SY.get_concrete_input(abstract_system, abstract_input)
            push!(concrete_input_traj, concrete_input)
        end
    end
    return ST.Control_trajectory(ST.Trajectory(concrete_state_traj), ST.Trajectory(concrete_input_traj))
end