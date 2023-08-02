# Defaut values for the ZMP based controller 

#-----------------------Global parameters-----------------------
zmax = 0.4225;                      # Initial height of the robot
Δz = 0.1 * zmax                     # Get down for about 10% of his initial height 
zc = 0.1924541398385977;            # Height of the CoM z-plane computed via URDF file and RigidBodyDynamics when the robot stands up   
ΔCoMz = 0.21916919708029195 - zc;   # Variation of the height computed via URDF file and RigidBodyDynamics    
Ts = 0.02;                          # Sampling frequency 
g = 9.81;                           # Gravity constant

#-----------------------Foot pattern parameters-----------------------
Lmax = 0.1;                         # Half of the step length 
θ_max = 12 * pi / 180;              # Max orientation change on the foot 
θ_0 = 0.0;                          # Initial orientation of the robot (robot frame)    
isLeftSupport = true;               # Initial support foot 

#-----------------------Path to follow definition-----------------------
# Straight path for 2D mode (default path) 
t = vec(0:100);                   # Parameter on the parametric equation of the path 
yPath = 1.18 .+ 0.0 .* t;
xPath = 0.01 * t;

initial_position = [xPath[1]; yPath[1]; θ_0]

#-----------------------ZMP generator block-----------------------
Tstep = 1;                  # Step period of the walking process
δ = 0.2;                    # Ratio between double support phase and the step period
Tdelay = 5;                 # Delay to allows the robot to reach his final height 
Twait = 5                   # Delay before walking 

#-----------------------Preview control parameter-----------------------
previewTime = 1             #  [s]
q_e = 1;                    # Weight of the error term in objective function 
r_u = 1e-6;                 # Weight of the input term in  objective function

#-----------------------ZMP preview control block-----------------------
xinit = [xPath[1]; 0.0; 0.0];    # Initial state on the X component of the ZMP 
yinit = [yPath[1]; 0.0; 0.0];    # Initial state on the Y component of the ZMP 

#-----------------------Swing Foot Trajectory Generator-----------------------
Tver = 0.0 * Tstep;     # Vertical time at the end of the step
hstep = 0.01;           # Maximal height of the foot position w.r.t world
