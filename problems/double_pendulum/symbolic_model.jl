using Symbolics

# Define symbolic variables for state, control input, and parameters
Symbolics.@variables x1 x2 x3 x4 u1 m1 m2 l1 l2 g
x = [x1, x2, x3, x4] # state vector
u = [u1]             # control vector

# Auxiliary variables
Δθ = x1 - x2
M = m1 + m2
α = m1 + m2 * sin(Δθ)^2

# Define the system dynamics
F1 = x3
F2 = x4
F3 =
    -sin(Δθ) * (m2 * l1 * x3^2 * cos(Δθ) + m2 * l2 * x4^2) -
    g * (M * sin(x1) - m2 * sin(x2) * cos(Δθ)) / (l1 * α) + u1
F4 =
    sin(Δθ) * (M * l1 * x3^2 + m2 * l2 * x4^2 * cos(Δθ)) +
    g * (M * sin(x1) * cos(Δθ) - M * sin(x2)) / (l2 * α)

# Define the vector field F
F = [F1, F2, F3, F4]

# Compute the Jacobian with respect to x
J = Symbolics.jacobian(F, x)

# Simplify the Jacobian if necessary (optional)
J_simplified = Symbolics.simplify.(J)

# Display the Jacobian
println("Jacobian of F with respect to x:")
display(J_simplified)

# Substitute specific values for parameters
params = Dict(m1 => 1.0, m2 => 1.0, g => 9.81, l1 => 1.0, l2 => 1.0)
J_substituted = [substitute(element, params) for element in J]

# Simplify the substituted Jacobian
J_simplified = Symbolics.simplify.(J_substituted)

# Display the simplified Jacobian
println("Simplified Jacobian with substituted parameters:")
display(J_simplified)
