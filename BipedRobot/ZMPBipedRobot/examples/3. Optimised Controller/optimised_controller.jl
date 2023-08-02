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
using DelimitedFiles

## Include and import the ZMP based controller
include(joinpath(@__DIR__, "..", "..", "src", "ZMPBipedRobot.jl"))
import .ZMPBipedRobot
const ZMProbot = ZMPBipedRobot

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = true; # Plot and save results
MODEL_2D = false;

SAVE_CSV = true;

local_dir = joinpath(@__DIR__, "..", "../")
saveFolder = local_dir * "examples/3. Optimised Controller"

# define file name to save
ref_fileName = saveFolder * "/walkingPattern_ref.csv"

###########################################################
#                    Simulation parameters                #
###########################################################
if MODEL_2D
    ## Straight path for 2D Robot Model
    t = vec(0:100)
    yPath = 1.18 .+ 0.0 .* t
    xPath = 0.01 * t
    θ_0 = 0.0      # initial orientation of the robot w.r.t the x-axis
    robot_model = "ZMP_2DBipedRobot.urdf"
else
    ## Circle path for 3D Robot Model
    t = vec(100:-1:75)
    xPath = -0 .- 1.18 * sin.(2 * pi / 100 .* t)
    yPath = 0 .+ 1.18 * cos.(2 * pi / 100 .* t)
    θ_0 = 0.0      # initial orientation of the robot w.r.t the x-axis
    robot_model = "ZMP_3DBipedRobot.urdf"
end
###########################################################
#                    ZMP based controller                 #
###########################################################
# Get the optimised paramters for a controller samplled at 50 Hz
candidates = readdlm(local_dir * "docs/2. Optimisation process/solutions_Fs50.txt", ',');

# Uncomment/comment the desired set of parameters 
best_key = 245 # Fast Trajectory
# best_key = 96 # Slow Trajectory

x = candidates[best_key, :]
Δz = x[1];
Tstep = x[2]
Lmax = x[3];
δ = x[4];
Tver = x[5];
hstep = x[6];

# Construct the best candidate
wo = ZMProbot.WalkingOptimization(Δz, Tstep, Lmax, δ, Tver, hstep, 0.4275, 5, NaN, NaN)
ZMProbot.computeAutoDefineParameters!(wo);
br = ZMProbot.defineBipedRobot(wo);

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

###########################################################
#                       Save Results                      #
###########################################################
if SAVE_CSV
    # Store into more convienant variables
    qr = ik.q_r
    ql = ik.q_l

    qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]
    tplot = reduce(vcat, zt.timeVec)

    CoM = reduce(hcat, ct.CoM)
    ZMP = reduce(hcat, zt.ZMP)

    # define headers for each column
    ref_header = ["time" "q1_l" "q1_r" "q2_l" "q2_r" "ZMPx" "ZMPy" "CoMx" "CoMy" "CoMz"]

    # combine arrays into a table
    ref_data = hcat(
        tplot,
        qref[:, 1],
        qref[:, 2],
        qref[:, 3],
        qref[:, 4],
        ZMP[1, :],
        ZMP[2, :],
        CoM[1, :],
        CoM[2, :],
        CoM[3, :],
    )
    ref_data2store = [ref_header; ref_data]

    # Save into a csv file
    open(ref_fileName, "w") do file
        CSV.write(ref_fileName, Tables.table(ref_data2store); delim = ',')
        return print("\nFile saved at $(ref_fileName)")
    end
end
