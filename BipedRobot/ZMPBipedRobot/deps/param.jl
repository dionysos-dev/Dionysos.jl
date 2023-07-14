#-----------------------Global parameters-----------------------
zmax = 0.4225;          # Initial height of the robot
Δz = 0.1 * zmax    
zc = 0.1924541398385977;         # Height of the CoM z-plane     
ΔCoMz = 0.21916919708029195 - zc;
Ts = 0.02;          # Sampling frequency 
g = 9.81;           # Gravity constant

#-----------------------Foot pattern parameters-----------------------
Lmax = 0.1;                # Half of the step length 
θ_max = 12* pi/180;        # Max orientation change on the foot 
d = 0.052 ;            # Centre-to-centre distance between the right foot and the left foot 
x0 = 0.0; y0 = 1.18; θ_0 = 0.0;   # Initial position of the robot (robot frame)    
initial_position = [x0; y0; θ_0]
isLeftSupport = true;      # Initial support foot 

#-----------------------Path to follow definition-----------------------
# Straight path
# pathType = "S"
t = vec(0 : 100);                # Parameter on the parametric equation of the path 
yPath = 1.18 .+ 0.0 .* t; 
xPath = 0.01 * t;

# Diagonal Straight path 
# pathType = "D"
# θ_0 = pi/4
# t = vec(0 : 100);              # Parameter on the parametric equation of the path 
# yPath = 0.01*t;
# xPath =  0.01*t;

# # Circle path
# pathType = "C"
# t = vec(100 : -1: 75); 
# xPath = -0 .- 1.18 * sin.(2*pi/100 .*t);
# yPath = 0 .+ 1.18 * cos.(2*pi/100 .*t);

# Straight + sinusoidale path
# t = vec(100 : - 1: 0); 
# yPath = 0.1 .- 0.1 * cos.(4*pi/100 .* t);
# xPath = 0.01 .* (100 .- t);

#-----------------------ZMP generator block-----------------------
Tstep = 1;                  # Step period of the walking process
δ = 0.2;                    # Ratio between double and single support phase
Tdelay = 5;                 # Delay at the start of the simulation
Twait = 0
#-----------------------Preview control parameter-----------------------
previewTime = 1             #  [s]
q_e = 1;
r_u = 1e-6;  

#-----------------------ZMP preview control block-----------------------
xinit = [x0;0.0;0.0];    # Initial state on the X component of the ZMP 
yinit = [y0;0.0;0.0];    # Initial state on the Y component of the ZMP 

#-----------------------Swing Foot Trajectory Generator-----------------------
Tver = 0.0 * Tstep;     # Vertical time at the end of the step
hstep = 0.01;           # maximal height of the foot position w.r.t world

#-----------------------Inverse Finematics block-----------------------
## TODO : make a code to get theses values from the URDF file 
L1 = 0.20125            # Length of the thigh
L2 = 0.172              # Length of the leg
offset_hip_to_motor = 0.04025;
offset_ankle_to_foot = 0.009;

