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
import .ZMPBipedRobot as ZMProbot

###########################################################
#                    Code parameters                      #
###########################################################
PLOT_RESULT = false;
MODEL_2D = true;

SAVE_CSV = false; 

saveFolder = "examples/1. Default Controller/results"

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
###########################################################
#                    ZMP based controller                 #
###########################################################

# Construct the biped robot which store the geomtrical propreties and the path wanted 
br = ZMProbot.BipedRobot(readFile = true,
                         URDFfileName = robot_model, 
                         paramFileName = "param.jl",
                         saveFolder = saveFolder
                         )
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

###########################################################
#                       Save Results                      #
###########################################################
if SAVE_CSV
    # Store into more convienant variables 
    qr = ik.q_r; 
    ql = ik.q_l; 
    qref = [ql[:, 1] qr[:, 1] ql[:, 2] qr[:, 2]]
    tplot = reduce(vcat, zt.timeVec)

    CoM = reduce(hcat, ct.CoM)
    ZMP = reduce(hcat, zt.ZMP) 

    # define headers for each column 
    # ref_header = ["time[s]" "q1_l[rad]" "q1_r[rad]"  "q2_l[rad]"  "q2_r[rad]"]
    ref_header = ["time" "q1_l" "q1_r"  "q2_l"  "q2_r" "ZMPx" "ZMPy" "CoMx" "CoMy" "CoMz"]

    # combine arrays into a table
    ref_data = hcat(tplot, qref[:, 1],  qref[:, 2],  qref[:, 3],  qref[:, 4], ZMP[1, :], ZMP[2, :], CoM[1, :], CoM[2, :], CoM[3, :])
    ref_data2store = [ref_header; ref_data]

    # Save into a csv file 
    open(ref_fileName, "w") do file
        CSV.write(ref_fileName,  Tables.table(ref_data2store), delim=',')
        print("\nFile saved at $(ref_fileName)")
    end
end 