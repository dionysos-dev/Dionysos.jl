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
include(joinpath(@__DIR__, "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot as ZMProbot

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = true;
SAVE_RESULT = false; 
ANIMATE_RESULT = true; 
MODEL_2D = false; 

saveFolder = "examples/3. Optimised Controller/results"

# define file name to save 
ref_fileName = saveFolder*"/walkingPattern_ref.csv"

###########################################################
#                    Simulation parameters                #
###########################################################
if MODEL_2D
    ## Straight path for 2D Robot Model 
    t = vec(0 : 100);              
    yPath = 1.18 .+ 0.0 .* t; 
    xPath = 0.01 * t;
    robot_model = "ZMP_2DBipedRobot.urdf"
else 
    ## Circle path for 3D Robot Model 
    t = vec(100 : -1: 75); 
    xPath = -0 .- 1.18 * sin.(2*pi/100 .*t);
    yPath = 0 .+ 1.18 * cos.(2*pi/100 .*t);
    robot_model = "ZMP_3DBipedRobot.urdf"
end 

# Plot Format 
lw = 2          # Plot line width 
dpi = 600       # dpi for saved figures
msw = 0.1       # Plot dot outerline width
ms = 2          # Plot dot size 

Δt = 1e-3       # Simulation step 

# Position control parameters
Kp = 10000.
Ki = 100.
Kd = 100.

# true : PD with dynamics compensation, false  : PD control 
ctrl = true 

###########################################################
#                    ZMP based controller                 #
###########################################################
# # Get the optimised paramters 
# candidates =  readdlm("examples/2. Optimisation process/results/solutions.txt", ',');
# best_key = 261

# Get the optimised paramters for a controller samplled at 50 Hz 
candidates =  readdlm("examples/2. Optimisation process/results/solutions_Fs50.txt", ',');
best_key = 247

x = candidates[best_key, :] 
# x[2]  = x[2] + 0.2
Δz = x[1] ;
Tstep = x[2];
Lmax = x[3];
δ = x[4];
Tver = x[5];
hstep = x[6];

# Construct the best candidate 
wo = ZMProbot.WalkingOptimization(  Δz, Tstep, Lmax, δ, Tver, hstep,
                                    0.4275, 5,
                                    NaN, NaN )
ZMProbot.computeAutoDefineParameters!(wo);
br = ZMProbot.defineBipedRobot(wo);               
br.saveFolder = saveFolder; 
br.xPath = xPath; br.yPath = yPath; 

# Construct the Preview Controller
pc = ZMProbot.PreviewController(br = br, check = PLOT_RESULT)

# Run the Foot Planer Algorithm and get the foot position 
fp = ZMProbot.FootPlanner(br = br, check = PLOT_RESULT)

# Get the ZMP reference trajectory 
zt = ZMProbot.ZMPTrajectory(br = br, fp = fp, check = PLOT_RESULT)

# Convert the ZMP reference trajectory into CoM trajectory
ct = ZMProbot.CoMTrajectory(br = br, pc = pc, zt = zt, check = PLOT_RESULT)

# Get the Swing Foot trajectory 
sf = ZMProbot.SwingFootTrajectory(br = br, fp = fp, zt = zt, check = PLOT_RESULT)

# Get the joint trajectory from the all path 
ik = ZMProbot.InverseKinematics(br = br, fp = fp, ct = ct, sf = sf, check = PLOT_RESULT)

# Store into more convienant variables 
qr = ik.q_r; 
ql = ik.q_l; 
qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]

CoM = reduce(hcat, ct.CoM)
ZMP = reduce(hcat, zt.ZMP)
tplot = reduce(vcat, zt.timeVec)

###########################################################
#                         Testbench                       #
###########################################################
stepL = reduce(hcat, sf.stepL)
stepR = reduce(hcat, sf.stepR)
plot(stepR[1, :], stepR[2, :], stepR[3, :]) 