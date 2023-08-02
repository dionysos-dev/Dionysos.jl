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
using DelimitedFiles

## Include and import the ZMP based controller 
include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot

const ZMProbot = ZMPBipedRobot

###########################################################
#                      Code parameters                    #
###########################################################
PLOT_RESULT = true;

ANIMATE_RESULT = true;

MODEL_2D = true;

# Folder to save figures 
local_dir = joinpath(@__DIR__, "..", "../")
saveFolder = local_dir * "docs/1. Default Controller"

###########################################################
#                    Simulation parameters                #
###########################################################

if MODEL_2D
    ## Straight path for 2D Robot Model 
    t = vec(0:100)
    yPath = 1.18 .+ 0.0 .* t
    xPath = 0.01 * t
    θ_0 = 0
    robot_model = "ZMP_2DBipedRobot.urdf"
else
    ## Circle path for 3D Robot Model 
    t = vec(100:-1:75)
    xPath = -0 .- 1.18 * sin.(2 * pi / 100 .* t)
    yPath = 0 .+ 1.18 * cos.(2 * pi / 100 .* t)
    θ_0 = 0
    robot_model = "ZMP_3DBipedRobot.urdf"
end

lw = 2          # Plot line width 
dpi = 600
msw = 0.1
ms = 2

# Simulation parameters
Δt = 1e-3       # Simulation step 

# Position control parameters
Kp = 10000.0
Ki = 100.0
Kd = 100.0

# true : PD with dynamics compensation, false  : random torque 
ctrl = true

###########################################################
#                    ZMP based controller                 #
###########################################################
# Construct the biped robot which store the geomtrical propreties and the path wanted 
br = ZMProbot.BipedRobot(;
    readFile = true,
    URDFfileName = robot_model,
    paramFileName = "param.jl",
)
br.xPath = xPath;
br.yPath = yPath;
br.initial_position = [xPath[1], yPath[1], θ_0]

# Construct the Preview Controller
pc = ZMProbot.PreviewController(; br = br, check = PLOT_RESULT)

# Run the Foot Planer Algorithm and get the foot position 
fp = ZMProbot.FootPlanner(; br = br, check = PLOT_RESULT)

# Get the ZMP reference trajectory 
zt = ZMProbot.ZMPTrajectory(; br = br, fp = fp, check = PLOT_RESULT)

# Convert the ZMP reference trajectory into CoM trajectory
ct = ZMProbot.CoMTrajectory(; br = br, pc = pc, zt = zt, check = PLOT_RESULT)

# Get the Swing Foot trajectory 
sf = ZMProbot.SwingFootTrajectory(; br = br, fp = fp, zt = zt, check = PLOT_RESULT)

# Get the joint trajectory from the all path 
ik = ZMProbot.InverseKinematics(; br = br, fp = fp, ct = ct, sf = sf, check = PLOT_RESULT)

# Store into more convienant variables 
qr = ik.q_r;
ql = ik.q_l;
qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]

CoM = reduce(hcat, ct.CoM)
ZMP = reduce(hcat, zt.ZMP)
tplot = reduce(vcat, zt.timeVec)

###########################################################
#                  Simulation environement                #
###########################################################
tend = tplot[end]       # Simulation time 
# Construct the robot in the simulation engine 
rs = ZMProbot.RobotSimulator(;
    fileName = robot_model,
    symbolic = false,
    add_contact_points = true,
    add_gravity = true,
    add_flat_ground = true,
);
# Generate the visualiser 
vis = ZMProbot.set_visulalizer(; mechanism = rs.mechanism)

# Intiial configuration 
boom = [0, 0]
actuators = [0, 0, 0, 0]
foot = [0, 0]
ZMProbot.set_nominal!(rs, vis, boom, actuators, foot)

# Simulate the robot 
controller! = ZMProbot.trajectory_controller!(rs, tplot, qref, Δt, Kp, Ki, Kd, ctrl)
ts, qs, vs = RigidBodyDynamics.simulate(rs.state, tend; Δt = Δt, controller!);

# Open the visulaiser and run the animation 
if ANIMATE_RESULT
    open(vis)
    MeshCatMechanisms.animate(vis, ts, qs)
end
###########################################################
#                      Plot results                       #
###########################################################

qsim = reduce(hcat, qs)';
vsim = reduce(hcat, vs)';
tsim = reduce(vcat, ts);
torque_sim = reduce(hcat, rs.torques)
CoMsim = reduce(hcat, rs.CoM)
dof_offset = length(boom);
ZMPsim = reduce(hcat, rs.ZMP)

if (length(torque_sim[1, :]) == length(tsim) + 1)
    len_sim = length(tsim)
    len_t = len_sim
elseif (length(torque_sim[1, :]) == length(tsim) - 1)
    len_t = length(torque_sim[1, :])
    len_sim = len_t
else
    len_t = length(tsim)
    len_sim = len_t
end

plt_θ = plot(;
    xlims = (0, tend),
    xlabel = L"$t$ [s]",
    legend = true,
    legendcolumns = 2,
    dpi = dpi,
    layout = (2, 2),
)
plt_ω = plot(;
    xlims = (0, tend),
    xlabel = L"$t$ [s]",
    ylabel = L"$\omega$ [rad/s]",
    layout = (2, 2),
    dpi = dpi,
)
plt_τ = plot(;
    xlims = (0, tend),
    dpi = dpi,
    xlabel = L"$t$ [s]",
    ylabel = L"$\tau$ [Nm]",
    layout = (2, 2),
)

for (side_idx, side) in enumerate(["Left", "Right"])
    for (joint, name) in enumerate(["Leg", "Knee"])
        plot!(
            plt_θ[joint, side_idx],
            tplot,
            qref[:, (2 * joint - 1) + (side_idx - 1)];
            lw = lw,
            label = "Reference",
            title = "$(side) $(name)",
        )
        plot!(
            plt_θ[joint, side_idx],
            tsim,
            qsim[:, dof_offset + (2 * joint - 1) + (side_idx - 1)];
            lw = lw,
            label = "Simulation",
        )
        str = latexstring("q_$(dof_offset + (2*joint - 1) + (side_idx - 1))") * " [rad]"
        ylabel!(plt_θ[joint, side_idx], str)

        str = latexstring("q̇_$(dof_offset + (2*joint - 1) + (side_idx - 1))") * " [rad/s]"
        ylabel!(plt_ω[joint, side_idx], str)

        str = latexstring("τ_$(dof_offset + (2*joint - 1) + (side_idx - 1))") * " [Nm]"
        ylabel!(plt_τ[joint, side_idx], str)
        plot!(
            plt_τ[joint, side_idx],
            tsim[1:len_t],
            torque_sim[dof_offset + (2 * joint - 1) + (side_idx - 1), 1:len_sim];
            lw = lw,
            label = false,
            title = "$(side) $(name)",
        )
        plot!(
            plt_ω[joint, side_idx],
            tsim,
            vsim[:, dof_offset + (2 * joint - 1) + (side_idx - 1)];
            lw = lw,
            label = false,
            title = "$(side) $(name)",
        )
        if (joint == 1)
            plot!(plt_θ[joint, side_idx]; legend = false)
        end
    end
end

plt_com = plot(; xlabel = L"$t$ [s]", xlims = (0, tend), layout = (2, 1), dpi = dpi)

plot!(
    plt_com[1],
    tsim[1:len_t],
    CoMsim[1, 1:len_sim];
    lw = lw,
    label = "Simulated",
    ylabel = L"$x$ [m]",
    title = "CoMx",
)
plot!(plt_com[1], tplot, CoM[1, :]; label = "Predicted", lw = lw)
plot!(
    plt_com[2],
    tsim[1:len_t],
    CoMsim[3, 1:len_sim];
    lw = lw,
    label = "Simulated",
    ylabel = L"$z$ [m]",
    title = "CoMz",
)
plot!(plt_com[2], tplot, CoM[3, :]; label = "Predicted", lw = lw)

ZMPref = reduce(hcat, zt.ZMP)
right_plot = reduce(hcat, fp.right)
left_plot = reduce(hcat, fp.left)
plt_ZMP = plot(;
    title = "ZMP in XY plane",
    xlabel = L"$x$ [m]",
    ylabel = L"$y$ [m]",
    legend = true,
    dpi = dpi,
)
scatter!(ZMPsim[1, :], ZMPsim[2, :]; ms = ms, msw = msw, label = "Measured ZMP")
plot!(ZMPref[1, :], ZMPref[2, :]; label = "Reference", lw = lw)
scatter!(left_plot[1, :], left_plot[2, :]; shape = :rect, label = "Left")
scatter!(right_plot[1, :], right_plot[2, :]; shape = :rect, label = "Right")

if PLOT_RESULT
    display(plt_θ)
    display(plt_τ)
    display(plt_ω)
    display(plt_com)
    display(plt_ZMP)
end
