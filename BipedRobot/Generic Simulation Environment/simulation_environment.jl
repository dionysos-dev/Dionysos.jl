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
import .GenericSimulationEnvironment as GSE

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = true;

ANIMATE_RESULT = true;
GRAVITY = true;
CONTACTS = false;
GROUND = true;

local_dir = joinpath(@__DIR__, "../../")
saveFolder = local_dir * "docs/6. Generic Simulation Environment"

refFolder = local_dir * "examples/6. Generic Simulation Environment/deps"

# define file name to open
ref_fileName = refFolder * "/walkingPattern_ref.csv"

###########################################################
#                    Simulation parameters                #
###########################################################
# Plot parameters
lw = 2
dpi = 600
msw = 0.1
ms = 2

# File Nme of the urdf file to open 
# fileName = "ZMP_3DBipedRobot.urdf"
fileName = "ZMP_2DBipedRobot_noContacts.urdf"

# Simulation parameters 
Δt = 1e-3       # Simulation step [s]
tend = 15        # Simulation time [s]

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

# open a reference joint trajectory 
csvpath() = joinpath(@__DIR__, "deps/", "$(ref_fileName)")
# Read the CSV file into a DataFrame
data = CSV.read(csvpath(), DataFrame; header = [2], delim = ',')
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
    hunt_crossley_hertz(α = α, k = k_n),
    ViscoelasticCoulombModel(µ, k_t, λ_t)
)

contacts_points = Array[]
push!(contacts_points, cp)
push!(contacts_points, cp1)
push!(contacts_points, cp2)
push!(contacts_points, cp3)
push!(contacts_points, cp4)

vr = GSE.VirtualRobot(
    fileName = fileName, 
    contactmodel =  contactmodel, 
    contacts_points = contacts_points, 
    add_contact_points = CONTACTS,
    add_flat_ground = GROUND,
    add_gravity =  GRAVITY
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
    flag = 1; 
    mechanism = vr.mechanism
    state = vr.state
    i =  0.0
    idx_sim = 0 
    
    # Do not touch the arguments of this function as it will be called at each sample 
    function controller!(τ, t, state)
        # Update the next torque as random value between 0 and 1 
        # rand!(τ)        
        # τ .= (τ .- 0.5)

        # Example of a basic controller 
        if t <= 1
            τ .= 0
            τ[3:3] .= 0 
        else 
            if isapprox(t, i * Ts, atol = Δt/10) 
                if (configuration(state)[3] <= -30 * pi/180)
                    τ[3:3] .= τ[3] + 0.002
                    flag = 1
                elseif  (configuration(state)[3] >= 30 * pi/180)
                    τ[3:3]  .= τ[3] - 0.002
                    flag = 0
                else 
                    if flag == 1 
                        τ[3:3] .= τ[3] + 0.002
                    else
                        τ[3:3] .= τ[3] - 0.002
                    end 
                end  
            end 
            if abs(τ[3]) >= 0.08
                τ[3:3] .= sign(τ[3]) *  0.08
            end     
        end 
        if isapprox(t, i * Ts, atol = Δt/10) 
            i = i + 1 
        end 
        if isapprox(t, idx_sim * Δt, atol = Δt/10)
            idx_sim = idx_sim + 1; 
            push!(vr.torques, copy(τ))
        end 
        return nothing
    end 
end 

###########################################################
#                  Simulation environement                #
###########################################################
GSE.set_initialbody!(vr)
# GSE.set_nominal!(vr, [0,0,30*pi/180,0,0,0])

# Define the controller 
controller! = define_controller!(vr, Δt, 0.2)

# Simulate the controller
ts, qs, vs = RigidBodyDynamics.simulate(vr.state, tend, Δt = Δt, controller!);

if ANIMATE_RESULT
    open(vr.vis) 
    MeshCatMechanisms.animate(vr.vis, ts, qs)
end 

###########################################################
#                      Plot results                       #
###########################################################
if PLOT_RESULT
    qsim = reduce(hcat, qs)';
    vsim = reduce(hcat, vs)';
    tsim = reduce(vcat, ts);

    torques_sim = reduce(hcat, vr.torques)'

    for q in range(1, length(vr.state.v))
        plt = plot(
            xlims = (0, tend),
            xlabel = L"$t$ [s]",
            legend = false,
            layout = (3, 1),
        )
        plot!(plt[1, 1], tsim,  qsim[:, q], ylabel = "q_$(q)")
        plot!(plt[2, 1], tsim,  vsim[:, q], ylabel = "̇q_$(q)")
        plot!(plt[3, 1], tsim,  torques_sim[:, q], ylabel = "τ_$(q)")
        
        display(plt)
    end 
end 