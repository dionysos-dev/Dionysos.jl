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
using DelimitedFiles,DataFrames

## Include and import the ZMP based controller 
include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = true;
SAVE_RESULT = false; 

ANIMATE_RESULT = true; 
GRAVITY = true; 
CONTACTS = false; 
GROUND = true; 

local_dir = joinpath(@__DIR__, "../../")
saveFolder = local_dir*"examples/4. Simulation Environment/results"
# refFolder = local_dir* "examples/1. Default Controller/results"
refFolder = local_dir*"examples/3. Optimised Controller/results"

# define file name to open
ref_fileName = refFolder*"/walkingPattern_ref_slow.csv"
# ref_fileName = refFolder*"/walkingPattern_ref_fast.csv"


if CONTACTS
    robot_model = "ZMP_3DBipedRobot.urdf"
else 
    robot_model = "ZMP_2DBipedRobot_noContacts.urdf"
end 
###########################################################
#                    Simulation parameters                #
###########################################################
# Plot parameters
lw = 2         
dpi = 600       
msw = 0.1 
ms = 2 

# Simulation parametersx
Δt = 1e-3       # Simulation step 

# Position control parameters
Kp = 10000.
Ki = 0.
Kd = 100.

# true : PD with dynamics compensation, false  : PD control 
ctrl = true 

###########################################################
#                  Simulation environement                #
###########################################################
data = ZMProbot.openCSV(ref_fileName)
q1_r = data.q1_r
q1_l = data.q1_l 
q2_r = data.q2_r
q2_l = data.q2_l
tplot = data.time 

ZMPx = data.ZMPx
ZMPy = data.ZMPy
CoMx = data.CoMx
CoMy = data.CoMy
CoMz = data.CoMz

qref = [q1_l q1_r q2_l q2_r]
ZMP = [ZMPx ZMPy]'
CoM = [CoMx CoMy CoMz]'
tend = tplot[end]       # Simulation time 

# Construct the robot in the simulation engine 
rs = ZMProbot.RobotSimulator(   fileName = robot_model,
                                symbolic = false, 
                                add_contact_points = CONTACTS, 
                                add_gravity = GRAVITY,
                                add_flat_ground = GROUND);
# Generate the visualiser 
vis = ZMProbot.set_visulalizer(mechanism = rs.mechanism)

# Intiial configuration
if CONTACTS
    boom = [ 0, 0 ]
else 
    boom = [ ]
end  
actuators = [ q1_l[1], q1_r[1], q2_l[1], q2_r[1] ]
foot = [0, 0]
ZMProbot.set_nominal!(rs, vis, boom, actuators, foot)

# Simulate the robot 
controller! = ZMProbot.trajectory_controller!(rs, tplot, qref, Δt, Kp, Ki, Kd, ctrl)
ts, qs, vs = RigidBodyDynamics.simulate(rs.state, tend, Δt = Δt,  controller!);

# Open the visulaiser and run the animation 
if ANIMATE_RESULT 
    open(vis) 
    MeshCatMechanisms.animate(vis, ts, qs)
end 
###########################################################
#                      Plot results                       #
###########################################################
qsim = reduce(hcat, qs)';
vsim =  reduce(hcat, vs)';
tsim = reduce(vcat, ts);
torque_sim =  reduce(hcat, rs.torques)
CoMsim =  reduce(hcat, rs.CoM)
dof_offset = length(boom); 
ZMPsim = reduce(hcat, rs.ZMP)

if(length(torque_sim[1, :]) == length(tsim) + 1) 
    len_sim =  length(tsim) 
    len_t = len_sim
elseif (length(torque_sim[1, :]) == length(tsim) - 1) 
    len_t =  length(torque_sim[1, :]) 
    len_sim = len_t
else 
    len_t = length(tsim)
    len_sim = len_t
end 

plt_θ = plot(   #title = "Joints angles",
                xlims = (0, tend), 
                xlabel = L"$t$ [s]",
                legend=true, legendcolumns=2, 
                #titlefont=font(fs, ff), 
                dpi = dpi,
                layout = (2,2))
plt_ω = plot(   #title = "Joints Velocity",
                xlims = (0, tend), 
                #ylims = (-pi, pi),
                xlabel = L"$t$ [s]", ylabel = L"$\omega$ [rad/s]", 
                layout = (2,2),
                dpi = dpi, )
plt_τ = plot(   #title = "Joints Torques",
                xlims = (0, tend), 
                #ylims = (-pi, pi),
                dpi = dpi, 
                xlabel = L"$t$ [s]", ylabel = L"$\tau$ [Nm]", 
                layout = (2,2))


for (side_idx, side) in enumerate(["Left", "Right"])
    for (joint, name) in enumerate(["Leg", "Knee"])
        plot!(plt_θ[joint, side_idx], tplot, qref[:,  (2*joint - 1) + (side_idx - 1)], lw = lw, label = "Reference", title = "$(side) $(name)")
        plot!(plt_θ[joint, side_idx], tsim, qsim[:,  dof_offset + (2*joint - 1) + (side_idx - 1)], lw = lw, label = "Simulation")
        str = latexstring("q_$(dof_offset + (2*joint - 1) + (side_idx - 1))" ) * " [rad/s]"
        ylabel!(plt_θ[joint, side_idx], str )

        str = latexstring("q̇_$(dof_offset + (2*joint - 1) + (side_idx - 1))" ) * " [rad/s]"
        ylabel!(plt_ω[joint, side_idx], str )

        str = latexstring("τ_$(dof_offset + (2*joint - 1) + (side_idx - 1))" ) * " [rad/s]"
        ylabel!(plt_τ[joint, side_idx], str )
        plot!(plt_τ[joint, side_idx], tsim[1:len_t], torque_sim[dof_offset + (2*joint - 1) + (side_idx - 1), 1:len_sim], lw = lw, label = false, title = "$(side) $(name)")
        plot!(plt_ω[joint, side_idx], tsim, vsim[:,   dof_offset + (2*joint - 1) + (side_idx - 1)], lw = lw, label = false, title = "$(side) $(name)")
        if (joint == 1)
            plot!(plt_θ[joint, side_idx], legend = false)
        end 
    end 
end 

plt_com = plot( xlabel = L"$t$ [s]",
            xlims = (0, tend), 
            layout = (2,1), dpi = dpi
            )

plot!(plt_com[1], tsim[1:len_t],  CoMsim[1, 1:len_sim], lw = lw, label = "Simulated", ylabel = L"$x$ [m]", title = "CoMx" )
plot!(plt_com[1], tplot, CoM[1, :], label = "Predicted",  lw = lw)
plot!(plt_com[2],tsim[1:len_t], CoMsim[3, 1:len_sim], lw = lw, label = "Simulated", ylabel = L"$z$ [m]", title = "CoMz")
plot!(plt_com[2], tplot, CoM[3, :], label = "Predicted",  lw = lw)

plt_ZMP = plot( title = "ZMP in XY plane",
                xlabel = L"$x$ [m]", ylabel = L"$y$ [m]",
                legend = true, dpi = dpi, 
                #aspect_ratio =:equal,
                #xlim = (-0.3, 1.5), ylim = (-0.3, 1.5)
                )
scatter!(ZMPsim[1, :], ZMPsim[2, :],  ms = ms, msw = msw, label = "Measured ZMP")
plot!(ZMP[1,:], ZMP[2,:] , label = "Reference",lw = lw)

if PLOT_RESULT 
    display(plt_θ) 
    display(plt_τ) 
    display(plt_ω) 
    display(plt_com)
    display(plt_ZMP) 
end 

if SAVE_RESULT
    savefig(plt_θ, saveFolder*"/sim_joints_angles.png")
    savefig(plt_ω, saveFolder*"/sim_joints_velocities.png")
    savefig(plt_τ, saveFolder*"/sim_joints_torques.png")
    savefig(plt_com,saveFolder*"/sim_CoMxz.png")
    savefig(plt_ZMP,saveFolder*"/sim_ZMPxy.png")
end 