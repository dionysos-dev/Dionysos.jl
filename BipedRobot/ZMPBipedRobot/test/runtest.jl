
# Test the ZMP based controller define in the param_test.jl file 
include("./ZMPBasedController/test_BipedRobot.jl")
include("./ZMPBasedController/test_PreviewControl.jl")
include("./ZMPBasedController/test_FootPlanner.jl")
include("./ZMPBasedController/test_ZMPTrajectoryGenerator.jl")
include("./ZMPBasedController/test_CoMTrajectoryGenerator.jl")
include("./ZMPBasedController/test_DirectKinematics.jl")
include("./ZMPBasedController/test_InverseKinematics.jl")
include("./ZMPBasedController/test_SwingFootTrajectoryGenerator.jl")

# Test the simulator 
include("./Simulation Environment/test_RobotSimulator.jl")
