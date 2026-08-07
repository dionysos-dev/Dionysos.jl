using StaticArrays, JuMP, Plots
using Symbolics, MathOptSymbolicAD

using Spot

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const MP = DI.Mapping
const AB = DI.Optim.Abstraction;

function integrator_model()
    model = Model(Dionysos.Optimizer)
    @variable(model, -2.0 <= x[1:2] <= 2.0)
    @variable(model, -1.0 <= u[1:2] <= 1.0)
    @constraint(model, ∂(x[1]) == u[1])
    @constraint(model, ∂(x[2]) == u[2])

    set_attribute(model, "jacobian_bound", u -> @SMatrix zeros(2, 2))
    set_attribute(model, "approx_mode", AB.UniformGridAbstraction.GROWTH)
    set_attribute(model, "time_step", 0.3)
    set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.2, 0.2)))
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))
    set_attribute(model, "print_level", 0)
    return model, x
end;

start_region = UT.box(SVector(-1.7, -1.7), SVector(-1.6, -1.6))
x0 = SVector(-1.65, -1.65);

model, x = integrator_model()

roomA = UT.set_union([
    UT.box(SVector(-1.0, 1.0), SVector(-0.3, 1.7)),
    UT.box(SVector(1.0, 1.0), SVector(1.7, 1.7)),
])
roomB = UT.box(SVector(-1.5, -1.2), SVector(-0.6, -0.2))
roomC = UT.box(SVector(1.0, -1.8), SVector(1.5, -1.1))
wall = UT.set_union([
    UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5)),
    UT.box(SVector(1.3, -0.5), SVector(2.0, 0.5)),
]);

@constraint(model, roomA_visited, x in Label(roomA))
@constraint(model, roomB_visited, x in Label(roomB))
@constraint(model, roomC_visited, x in Label(roomC))
@constraint(model, wall_hit, x in Label(wall; semantics = MP.OUTER))

@constraint(model, x in Start(start_region));

@specification(
    model,
    ltl"G(!wall_hit) & F(roomA_visited & F(roomB_visited & F(roomC_visited)))"
)

optimize!(model);

termination_status(model)

sequencing_problem = get_attribute(model, "concrete_problem");

sequencing_traj = Dionysos.simulate(model, x0; nsteps = 80);

fig1 = plot(; aspect_ratio = :equal, title = "visit A, then B, then C")
plot!(
    sequencing_problem;
    ap_colors = Dict(
        :roomA_visited => :green,
        :roomB_visited => :cyan,
        :roomC_visited => :orange,
        :wall_hit => :black,
    ),
)
plot!(sequencing_traj; ms = 1.5, color = :blue)

const OPDS = DI.Optim.DiscreteSystems;

function until_step(q::Int, ap::Tuple{Vararg{Symbol}})
    barrier = :barrier_hit in ap
    first_stop = :first_stop_reached in ap
    second_stop = :second_stop_reached in ap
    danger = :danger_zone in ap

    barrier && return 0          ## the barrier is fatal at every point in the task
    q == 0 && return 0           ## dead absorbs
    q == 3 && return 3           ## done absorbs

    if q == 1
        # Before the first stop the danger zone is free to cross; arriving at the first
        # stop while already standing in the second finishes the task outright.
        first_stop || return 1
        return second_stop ? 3 : 2
    end

    # q == 2: the until is running — reach the second stop, and no danger before then.
    second_stop && return 3
    return danger ? 0 : 2
end

until_monitor = OPDS.FunctionMonitor(1, Set([3]), until_step);

model2, x2 = integrator_model()

first_stop = UT.box(SVector(1.0, 1.0), SVector(1.7, 1.7))
second_stop = UT.box(SVector(-1.5, -1.2), SVector(-0.6, -0.2))
barrier = UT.box(SVector(-1.8, 0.0), SVector(-0.6, 1.0))
danger = UT.set_union([
    UT.box(SVector(-0.5, -0.5), SVector(0.5, 0.5)),
    UT.box(SVector(1.3, -0.5), SVector(2.0, 0.5)),
]);

@constraint(model2, first_stop_reached, x2 in Label(first_stop))
@constraint(model2, second_stop_reached, x2 in Label(second_stop))
@constraint(model2, barrier_hit, x2 in Label(barrier; semantics = MP.OUTER))
@constraint(model2, danger_zone, x2 in Label(danger; semantics = MP.OUTER))

@constraint(model2, x2 in Start(start_region));

@specification(model2, until_monitor)

optimize!(model2);

termination_status(model2)

until_problem = get_attribute(model2, "concrete_problem")
until_traj = Dionysos.simulate(model2, x0; nsteps = 80);

fig2 = plot(; aspect_ratio = :equal, title = "reach A, then avoid danger until B")
plot!(
    until_problem;
    ap_colors = Dict(
        :first_stop_reached => :green,
        :second_stop_reached => :cyan,
        :barrier_hit => :black,
        :danger_zone => :red,
    ),
)
plot!(until_traj; ms = 1.5, color = :blue)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
