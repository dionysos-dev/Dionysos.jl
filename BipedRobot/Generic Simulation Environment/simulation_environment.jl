## Import useful packages 
using LinearAlgebra
using StaticArrays
using StructArrays
using Plots
using RigidBodyDynamics
using RigidBodyDynamics.Contact
using StaticArrays
using Symbolics
using MeshCat, MeshCatMechanisms, Blink
using MechanismGeometries
using LaTeXStrings
using Tables, CSV
using DelimitedFiles, DataFrames
using Random

## Include and import the ZMP based controller 
include(joinpath(@__DIR__, "utils/", "GenericSimulationEnvironment.jl"))
import .GenericSimulationEnvironment

const GSE = GenericSimulationEnvironment

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = true;

ANIMATE_RESULT = true;
GRAVITY = true;
CONTACTS = false;
GROUND = true;

local_dir = joinpath(@__DIR__)

refFolder =
    joinpath(local_dir, "..", "ZMPBipedRobot", "examples", "3. Optimised Controller")

# define file name to open
ref_fileName = joinpath(refFolder, "walkingPattern_ref.csv")

###########################################################
#                    Simulation parameters                #
###########################################################
# Plot parameters
lw = 2
dpi = 600
msw = 0.1
ms = 2

# File Nme of the urdf file to open 
fileName = "ZMP_2DBipedRobot_noContacts.urdf"

# Simulation parameters 
Δt = 1e-3       # Simulation step [s]
tend = 20        # Simulation time [s]

# Contact Points location in the ankle foot frame  
cp = [0.0, 0.0, -0.009]
cp1 = [0.035, 0.02, -0.009]
cp2 = [0.035, -0.02, -0.009]
cp3 = [-0.035, 0.02, -0.009]
cp4 = [-0.035, -0.02, -0.009]

# Normal Contact Model (Hunt-Crossley)
k_n = 50e3      # Stiffness constant [N/m]
α = 0.2         # 1st order coefficient of the linerised coefficient of restitution [s/m]

# Tangent Contact Model (Coulomb)
μ = 0.8         # Coefficient of friction between the foot and the ground [.] 
k_t = 20e3      # Stiffness constant [N/m]
λ_t = 100.0     # Damping constant [kg/s]

# Read the CSV file into a DataFrame]
data = CSV.read(ref_fileName, DataFrame; header = [2], delim = ',')
q1_r = data.q1_r
q1_l = data.q1_l
q2_r = data.q2_r
q2_l = data.q2_l
tref = data.time

qref = [q1_l q1_r q2_l q2_r]

###########################################################
#                     Robot Definition                    #
###########################################################

# Contact Model definition
contactmodel = SoftContactModel(
    hunt_crossley_hertz(; α = α, k = k_n),
    ViscoelasticCoulombModel(µ, k_t, λ_t),
)

contacts_points = Array[]
push!(contacts_points, cp)
push!(contacts_points, cp1)
push!(contacts_points, cp2)
push!(contacts_points, cp3)
push!(contacts_points, cp4)

vr = GSE.VirtualRobot(;
    fileName = fileName,
    contactmodel = contactmodel,
    contacts_points = contacts_points,
    add_contact_points = CONTACTS,
    add_flat_ground = GROUND,
    add_gravity = GRAVITY,
)

###########################################################
#                   Controller Definition                 #
###########################################################
# This controller is called 4 time for one sample Δt 

function define_controller!(
    vr::GSE.VirtualRobot,
    Δt,
    Ts,
    # Feel free to add more arguments to complexify the controller
)
    flag = 1
    mechanism = vr.mechanism
    state = vr.state
    i = 0.0
    idx_sim = 0

    # Do not touch the arguments of this function as it will be called at each sample 
    function controller!(τ, t, state)
        # Example of a basic controller for one joint 
        if t <= 1
            τ .= 0
        else
            if isapprox(t, i * Ts; atol = Δt / 10)
                if (configuration(state)[2] <= -30 * pi / 180)
                    τ[2:2] .= τ[2] + 0.02
                    flag = 1
                elseif (configuration(state)[2] >= 30 * pi / 180)
                    τ[2:2] .= τ[2] - 0.02
                    flag = 0
                else
                    if flag == 1
                        τ[2:2] .= τ[2] + 0.02
                    else
                        τ[2:2] .= τ[2] - 0.02
                    end
                end
            end
            # Saturation 
            if abs(τ[2]) >= 0.5
                τ[2:2] .= sign(τ[2]) * 0.5
            end
        end
        if isapprox(t, i * Ts; atol = Δt / 10)
            i = i + 1
        end
        if isapprox(t, idx_sim * Δt; atol = Δt / 10)
            idx_sim = idx_sim + 1
            push!(vr.torques, copy(τ))
        end
        return nothing
    end
end

###########################################################
#                  Simulation environement                #
###########################################################
GSE.set_initialbody!(vr)

# Define the controller 
controller! = define_controller!(vr, Δt, 0.2)

# Simulate the controller
ts, qs, vs = RigidBodyDynamics.simulate(vr.state, tend, controller!; Δt = Δt);

if ANIMATE_RESULT
    open(vr.vis)
    MeshCatMechanisms.animate(vr.vis, ts, qs)
end

###########################################################
#                      Plot results                       #
###########################################################
if PLOT_RESULT
    qsim = reduce(hcat, qs)'
    vsim = reduce(hcat, vs)'
    tsim = reduce(vcat, ts)

    torques_sim = reduce(hcat, vr.torques)'

    for q in 1:length(vr.state.v)
        p = plot(;
            xlims = (0, tend),
            xlabel = "\$$t [s]\$",
            legend = false,
            layout = (3, 1),
        )
        plot!(p[1, 1], tsim, qsim[:, q]; ylabel = "\$q_{$q}\$")
        plot!(p[2, 1], tsim, vsim[:, q]; ylabel = "\$q_{$q}\$")
        plot!(p[3, 1], tsim, torques_sim[:, q]; ylabel = "\$\\tau_{$q}\$")
        display(p)
    end
end
