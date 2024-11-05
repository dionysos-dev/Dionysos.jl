using Test
# # Example: Thermostat problem solved by [Uniform grid abstraction](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/src/manual/manual.md#solvers).
# 
# Let us consider the 1-dimensional state space control system of the form
# ```math
# \dot{x} = f(T, u)
# ```
# with $T=T(t)$ the temperature of the room at time $t$ and $u=u(t)$ the heating intensity, where $0 \leq u(t) \leq 1$. A value of 1 represents maximum heating, and 0 represents no heating.
# $f: \mathbb{R} × U \to \mathbb{R}$. The system can operate in two modes: heating (mode=1) and idle (mode=0).
# 1. In heating mode, the temperature increases at a rate proportional to the heating intensity (the control input $u(t)$).
# ```math
# f(T, u) = k \cdot u(t)
# ```
# 2. In idle mode, the temperature decreases (or remains constant if $\alpha=0$) in this mode due to the abient cooling.
# ```math
# f(T, u) = -\alpha
# ```
# where $k$ is Maximum heating rate (degrees per unit time).
#
# The control objective is to maintain the temperature of the room at a desired reference temperature $T_r$ with a tolerance of $\delta$. This introduces a constraint on the control input $u(t)$.
# ```math
# |T(t) - T_r| \leq \delta
# ``` 
# or equivalently
# ```math
# T_r - \delta \leq T(t) \leq T_r + \delta
# ```
#
# In this example, we will consider the hybrid system where the system can switch between the two modes.

# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl).
using StaticArrays, Plots

# At this point, we import Dionysos and JuMP.
using Dionysos, JuMP

# We define the parameters of the system.
k = 0.5 # Maximum heating rate (degrees per unit time)
α = 0.1 # Ambient cooling rate (degrees per unit time)
T_r = 20.0 # Reference temperature
δ = 1.0 # Tolerance

# ### Definition of the problem

# We define the problem using JuMP as follows.
# We first create a JuMP model:
model = Model(Dionysos.Optimizer)

# Define the state variables: T(t)
T_low, T_upp = 0.0, 30.0
T_start = 10.0
@variable(model, T_low <= T <= T_upp, start = T_start)

# Define the control variables: u(t)
@variable(model, 0 <= u <= 1)

# Define the mode variable: mode(t)
@variable(model,  0 <= mode <= 1, integer=true)

# Define the dynamics
#@constraint(model, ∂(T) == (1 - mode) * (-α) + mode * k * u)

@constraint(model, mode == 1 => {∂(T) == k * u} )
@constraint(model, mode == 0 => {∂(T) == -α} )

#@constraint(model, ifelse(mode == 1, ∂(T) == k * u, ∂(T) == -α))

@constraint(model, T_r - δ <= T)
@constraint(model, T <= T_r + δ)

# Define the initial set (from the other PR to make this work)
#@constraint(model, start(T) in MOI.Interval(9.5, 10.5))

# Define the target set
@constraint(model, final(T) in MOI.Interval(19.0, 21.0))

# Define the time step
set_attribute(model, "time_step", 0.1)

# ### Definition of the abstraction

# Definition of the grid of the state-space on which the abstraction is based (origin `x0` and state-space discretization `h`):
x0 = SVector(0.0);
h = SVector(0.1);
set_attribute(model, "state_grid", Dionysos.Domain.GridFree(x0, h))

# Definition of the grid of the input-space on which the abstraction is based (origin `u0` and input-space discretization `h`):
u0 = SVector(0.0);
h = SVector(0.3);
set_attribute(model, "input_grid", Dionysos.Domain.GridFree(u0, h))

optimize!(model)
