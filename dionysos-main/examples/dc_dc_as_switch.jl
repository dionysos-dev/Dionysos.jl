using CoreControls
using Solver

# Create a Control model for the switching system
model = Control(
    name="SwitchingSystem",
    problem_type=Safety(),
    system_type=Hybrid()  # Hybrid system type due to mode switching
)

# Parameters
@parameter(model, "r_l", 1.0)
@parameter(model, "x_l", 1.0)
@parameter(model, "r_c", 1.0)
@parameter(model, "x_c", 1.0)
@parameter(model, "r_0", 1.0)
@parameter(model, "x_r", 1.0)
@parameter(model, "u_s", 1.0)
@parameter(model, "tau", 0.1)  # Sampling time
@parameter(model, "N", 10)  # Horizon

# State Variables (i_l and v_c)
@variable(model, "i_l", StateVar(), Reals())
@variable(model, "v_c", StateVar(), Reals())

# Control Input (Switching Signal as a mode variable)
@variable(model, "switch_signal", ModeVar(), Integers(), bounds=(1, 2))

# System Dynamics (Mode 1 and Mode 2)
A1 = [
    [-r_l / x_l  0; 
     0           -1 / x_c  + 1 / (r_0 + r_c)]
]
A2 = [
    [-1 / x_l  (r_l + r_0 * x_r / (r_0 + r_c)) / x_l;
     -1 / x_c  r_0 / (r_0 + r_c)]
]
b = [u_s / x_l; 0]

for t in 0:N-1
    @constraint(model, :(dot(i_l, t, tau, model) == ifelse(switch_signal == 1, A1[1,1]*i_l + A1[1,2]*v_c + b[1], A2[1,1]*i_l + A2[1,2]*v_c + b[1])))
    @constraint(model, :(dot(v_c, t, tau, model) == ifelse(switch_signal == 1, A1[2,1]*i_l + A1[2,2]*v_c + b[2], A2[2,1]*i_l + A2[2,2]*v_c + b[2])))
end

# Safety Constraint: Keep system within safety region
@constraint(model, :(i_l^2 + v_c^2 <= 1))  # Example safety constraint

# Objective: Minimize deviation from reference
@objective(model, Minimize(), :(abs(i_l - 1.0) + abs(v_c - 1.0)))

# Solve the problem
algorithm = UniformGridAlgorithm()
solution = model.solve(algorithm; horizon=N, Δt=tau)

# Visualize the solution (Optional)
# print(solution)
