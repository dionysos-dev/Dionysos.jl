using StaticArrays, Plots

using Dionysos, JuMP

model = Model(Dionysos.Optimizer);

x_low, x_upp = [0.0, 0.0, -pi - 0.4], [4.0, 10.0, pi + 0.4]
x_start = [0.4, 0.4, 0.0]
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = x_start[i]);

@variable(model, -1 <= u[1:2] <= 1);

@expression(model, α, atan(tan(u[2]) / 2))

@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3]) == u[1] * tan(u[2]))

x_target = [3.3, 0.5, 0]

@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
@constraint(model, final(x[2]) in MOI.Interval(0.3, 0.8))
@constraint(model, final(x[3]) in MOI.Interval(-100.0, 100.0))

x1_lb = [1.0, 2.2, 2.2]
x1_ub = [1.2, 2.4, 2.4]
x2_lb = [0.0, 0.0, 6.0]
x2_ub = [9.0, 5.0, 10.0];

for i in eachindex(x1_ub)
    @constraint(
        model,
        x[1:2] ∉ MOI.HyperRectangle([x1_lb[i], x2_lb[i]], [x1_ub[i], x2_ub[i]])
    )
end

function jacobian_bound(u)
    β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
    return StaticArrays.SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
end
set_attribute(model, "jacobian_bound", jacobian_bound)
set_attribute(model, "time_step", 0.3)
set_attribute(
    model,
    "approx_mode",
    Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
)
set_attribute(model, "efficient", true)

x0 = SVector(0.0, 0.0, 0.0);
hx = SVector(0.2, 0.2, 0.2);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, hx))

u0 = SVector(0.0, 0.0);
hu = SVector(0.3, 0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, hu))

optimize!(model);

abstract_system = get_attribute(model, "abstract_system");
abstract_problem = get_attribute(model, "abstract_problem");
abstract_controller = get_attribute(model, "abstract_controller");
concrete_controller = get_attribute(model, "concrete_controller")
concrete_problem = get_attribute(model, "concrete_problem");
abstract_value_function = get_attribute(model, "abstract_value_function")
concrete_system = concrete_problem.system
abstraction_time =
    MOI.get(model, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
println("Time to construct the abstraction: $(abstraction_time)")
abstract_problem_time =
    MOI.get(model, MOI.RawOptimizerAttribute("abstract_problem_time_sec"))
println("Time to solve the abstract problem: $(abstract_problem_time)")
total_time = MOI.get(model, MOI.RawOptimizerAttribute("solve_time_sec"))
println("Total time: $(total_time)")

nstep = 100
function reached(x)
    if x ∈ concrete_problem.target_set
        return true
    else
        return false
    end
end

x0 = SVector(0.4, 0.4, 0.0)
control_trajectory = Dionysos.System.get_closed_loop_trajectory(
    get_attribute(model, "discrete_time_system"),
    concrete_controller,
    x0,
    nstep;
    stopping = reached,
)

using Plots

fig = plot(; aspect_ratio = :equal);

plot!(concrete_system.X; color = :grey, opacity = 1.0, label = "");

plot!(abstract_system; value_function = abstract_value_function);

plot!(concrete_problem.initial_set; color = :green, opacity = 0.2, label = "Initial set");
plot!(
    concrete_problem.target_set;
    dims = [1, 2],
    color = :red,
    opacity = 0.5,
    label = "Target set",
);

plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.initial_set);
    color = :green,
);
plot!(
    Dionysos.Symbolic.get_domain_from_states(abstract_system, abstract_problem.target_set);
    color = :red,
);

plot!(control_trajectory; ms = 2.0, arrows = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
