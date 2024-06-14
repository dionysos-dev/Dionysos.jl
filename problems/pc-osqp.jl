module PCOSQP
using Test
# First, let us import [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl).

using StaticArrays
using MathematicalSystems
using LinearAlgebra
using BlockDiagonals

# At this point, we import the useful Dionysos sub-module for this problem:
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const SY = DI.Symbolic

# Define the structure for the ADP parameters:
struct ADPParameters
    gamma::Float64
    N_adp::UnitRange{Int}
    delta::Float64
    N_tail::Int

    function ADPParameters()
        return new(0.95, 1:5, 5.5, 50)
    end
end

# Define the structure for the simulation parameters:
struct SimulationParameters
    Ts::Float64
    freq::Float64
    torque::Float64
    t0::Float64
    init_periods::Int
    sim_periods::Int
    flag_steady_trans::Int

    function SimulationParameters()
        Ts = 25.0e-06            # Sampling time
        freq = 50.0               # Switching frequency
        torque = 1.0              # Desired torque
        t0 = 0.0                 # Initial time
        init_periods = 1         # Number of periods to settle before simulation
        sim_periods = 2          # Numer of simulated periods
        flag_steady_trans = 0    # Flag Steady State (0) or Transients (1)

        return new(Ts, freq, torque, t0, init_periods, sim_periods, flag_steady_trans)
    end
end

# Define the structure for the parameters:
struct Parameters
    wr::Float64
    is::Float64
    Rs::Float64
    Rr::Float64
    Xls::Float64
    Xlr::Float64
    Xm::Float64
    Vdc::Float64
    xc::Float64
    Vrat::Float64
    Irat::Float64
    RealPower::Float64
    ApparentPower::Float64
    RotationalSpeed::Float64
    iVdc::Float64
    k1::Float64
    k2::Float64
    fsw_des::Float64
    kT::Float64
    P::SMatrix
    invP::SMatrix
    adp_params::ADPParameters
    sim_params::SimulationParameters

    function Parameters()
        ## P is a matrix of 2 rows and 3 columns with the following values:
        Park = (2.0 / 3.0) * @SMatrix [1.0 -0.5 -0.5; 0.0 sqrt(3.0)/2.0 -sqrt(3.0)/2.0]

        invPark = @SMatrix [1.0 0.0; -0.5 sqrt(3.0)/2.0; -0.5 -sqrt(3.0)/2.0]

        return new(
            0.9911,  # wr - Nominal speed
            1.0, #is - pu
            0.0108, #Rs - pu
            0.0091, #Rr - pu
            0.1493, #Xls - pu
            0.1104, #Xlr - pu
            2.3489, #Xm - pu
            1.930, #Vdc - pu
            11.769, #xc - pu
            3300, #Vrat - Voltage - V
            356, #Irat - Current - A
            1.587, #RealPower - MW
            2.035, #ApparentPower - MVA
            596, #RotationalSpeed - rpm
            5.2 * 1000, #iVdc - kV
            #switching filter parameters
            0.8e03, #k1
            0.8e03, #k2
            300, #fsw_des: desired switching frequency
            # Torque constant to fix [pu] units
            1.2361, #kT
            # Define Park transformation and its inverse
            Park, # P
            invPark, # invP
            ADPParameters(), # ADP parameters
            SimulationParameters(), # Simulation parameters
        )
    end
end

# Define the structure for the derived parameters:
struct DerivedParameters
    Vb::Float64
    Ib::Float64
    Fb::Float64
    Xs::Float64
    Xr::Float64
    taus::Float64 # tau of the stator
    taur::Float64 # tau of the rotor
    Dval::Float64
    Tspu::Float64
    omegar::Float64
    F::SMatrix
    G_first::SMatrix
    A_phys::SMatrix
    B_phys::SMatrix
    A_osc::SMatrix
    B_osc::SMatrix
    A_prev::SMatrix
    B_prev::SMatrix
    a1::Float64
    a2::Float64
    A_sw::SMatrix
    B_sw::SMatrix
    A_fsw::SMatrix
    A::SMatrix
    B::SMatrix
    C::SMatrix
    W::SMatrix
    G::SMatrix
    T::SMatrix
    M::SMatrix

    function DerivedParameters(this::Parameters)
        # Definition of the parameters of the system:

        #- From [9].sect I
        Vb = sqrt(1.0 / 3.0) * this.Vrat
        Ib = sqrt(2) * this.Irat
        Fb = this.sim_params.freq
        #- From [9].sect II.B
        ## Inductances
        Xs = this.Xls + this.Xm
        Xr = this.Xlr + this.Xm
        Dval = 0.6266 # Xs * Xr - this.Xm^2

        ## Time constants
        Tau_s = Xr * Dval / (this.Rs * Xr^2 + this.Rr * this.Xm^2) # not sure about this
        Tau_r = Xr / this.Rr

        Tspu = this.sim_params.Ts * this.wr
        omega_r = this.wr

        # Physical system matrices
        F = @SMatrix [
            -1.0/Tau_s 0.0 this.Xm/(Tau_r * Dval) omega_r * this.Xm/Dval
            0.0 -1.0/Tau_s -omega_r * this.Xm/Dval this.Xm/(Tau_r * Dval)
            this.Xm/Tau_r 0.0 -1.0/Tau_r -omega_r
            0.0 this.Xm/Tau_r omega_r -1.0/Tau_r
        ]

        # Define the matrix G
        Gtmp = @SMatrix [1.0 0.0; 0.0 1; 0.0 0.0; 0.0 0.0]
        G_first = (Xr / Dval * this.Vdc / 2.0) * Gtmp * this.P

        ##using LinearAlgebra, StaticArrays

        # Discretize physical system
        A_phys = exp(F .* Tspu)
        B_phys = -(inv(F) * (I - A_phys) * G_first)

        # Concatenate oscillating states
        A_osc = @SMatrix [
            cos(Tspu) -sin(Tspu)
            sin(Tspu) cos(Tspu)
        ]
        B_osc = zeros(Float64, 2, 3)

        # Concatenate previous input as a state
        A_prev = zeros(Float64, 3, 3)
        B_prev = I(3)

        # Concatenate filter states
        a1 = 1.0 - 1.0 / this.k1
        a2 = 1.0 - 1.0 / this.k2

        A_sw = @SMatrix [a1 0.0; (1.0-a1) a2]
        B_sw =
            (1.0 / this.fsw_des * 1.0 / 12 * (1 - a1) / this.sim_params.Ts) *
            @SMatrix [1.0 1.0 1.0; 0.0 0.0 0.0]

        # Concatenate switching frequency state
        A_fsw = @SMatrix [1.0]

        # Generate complete system
        A = BlockDiagonal([A_phys, A_osc, A_prev, A_sw, A_fsw])
        B = [
            B_phys zeros(size(B_phys, 1), 3)
            B_osc zeros(2, 3)
            B_prev zeros(3, 3)
            zeros(2, 3) B_sw
            zeros(1, 6)
        ]

        C = @SMatrix [
            1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -this.adp_params.delta this.adp_params.delta
        ]

        # Extract previous input from data
        W = [zeros(Float64, 3, 6) I(3) zeros(Float64, 3, 3)]

        # Matrix to extract original input and auxiliary one from state
        G = [I(3) zeros(Float64, 3, 3)]
        T = [zeros(Float64, 3, 3) I(3)]

        # Generate matrix of possible input combinations
        M = zeros(Float64, 3, 27)
        M[1, :] = [fill(-1.0, 9); fill(0.0, 9); fill(1.0, 9)]
        M[2, :] = repeat([fill(-1.0, 3); fill(0.0, 3); fill(1.0, 3)]; outer = [1, 3])
        M[3, :] = repeat([-1.0, 0.0, 1.0]; inner = [1], outer = [9])

        return new(
            Vb,
            Ib,
            Fb,
            Xs,
            Xr,
            Tau_s,
            Tau_r,
            Dval,
            Tspu,
            omega_r,
            F,
            G_first,
            A_phys,
            B_phys,
            A_osc,
            B_osc,
            A_prev,
            B_prev,
            a1,
            a2,
            A_sw,
            B_sw,
            A_fsw,
            A,
            B,
            C,
            W,
            G,
            T,
            M,
        )

        #new(Vb, Ib, Fb, Xs, Xr, Tau_s, Tau_r, Dval, Tspu, omega_r, F, G_first, A_phys, B_phys, A_osc, B_osc, A_prev, B_prev, a1, a2, A_sw, B_sw, A_fsw, A, B, C, W, G, T, M)
    end
end

# Define the structure for the simulation results:
struct SimulationResults
    X::Any
    U::Any
    Y_phase::Any
    Y_star_phase::Any
    T_e::Any
    T_e_des::Any
    solve_times::Any

    function SimulationResults(X, U, Y_phase, Y_star_phase, T_e, T_e_des, solve_times)
        return new(X, U, Y_phase, Y_star_phase, T_e, T_e_des, solve_times)
    end
end

function Base.show(io::IO, emp::SimulationResults)
    println(io, "SimulationResults")
    println(io, "X: ", emp.X)
    println(io, "U: ", emp.U)
    println(io, "Y_phase: ", emp.Y_phase)
    println(io, "Y_star_phase: ", emp.Y_star_phase)
    println(io, "T_e: ", emp.T_e)
    println(io, "T_e_des: ", emp.T_e_des)
    return println(io, "solve_times: ", emp.solve_times)
end

test_parameters = Parameters()
test_derived_parameters = DerivedParameters(test_parameters)

function dynamicofsystem(
    params::Parameters = test_parameters,
    dparams::DerivedParameters = test_derived_parameters,
)

    # Get parameters
    taus = dparams.taus
    taur = params.taur
    D = params.Dval
    omegar = params.wr
    Vdc = params.Vdc
    Xm = params.Xm
    Xr = dparams.Xr
    P = params.P
    Tspu = dparams.Tspu
    k1 = params.k1
    k2 = params.k2
    Ts = params.sim_params.Ts

    # Define the system dynamics

    # Definition of the growth bound function of $f$:
    A2_abs = SMatrix{2, 2}(
        -(rL + r0 * rC / (r0 + rC)) / xL,
        5.0 * r0 / (r0 + rC) / xC,
        r0 / (r0 + rC) / xL / 5.0,
        -1.0 / xC / (r0 + rC),
    )

    L_growthbound = let A1 = A1, A2_abs = A2_abs
        u -> u[1] == 1 ? A1 : A2_abs
    end

    DF_sys = let A1 = A1, A2 = A2
        (x, u) -> u[1] == 1 ? A1 : A2
    end
    bound_DF(u) = 0.0
    bound_DDF(u) = 0.0
    return F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF
end

function system(;
    sysnoise = SVector(0.0, 0.0),
    measnoise = SVector(0.0, 0.0),
    tstep = 0.5,
    nsys = 5,
    _X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85)),
    _U_ = UT.HyperRectangle(SVector(1), SVector(2)),
    xdim = 2,
    udim = 1,
    approx_mode = "growth",
)

    # Definition of the dynamics functions $f_p$ of the system:
    F_sys, L_growthbound, ngrowthbound, DF_sys, bound_DF, bound_DDF = dynamicofsystem()
    contsys = nothing
    if approx_mode == "growth"
        contsys = ST.NewControlSystemGrowthRK4(
            tstep,
            F_sys,
            L_growthbound,
            sysnoise,
            measnoise,
            nsys,
            ngrowthbound,
        )
    elseif approx_mode == "linearized"
        contsys = ST.NewControlSystemLinearizedRK4(
            tstep,
            F_sys,
            DF_sys,
            bound_DF,
            bound_DDF,
            measnoise,
            nsys,
        )
    end
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        contsys,
        xdim,
        udim,
        _X_,
        _U_,
    )
end

function problem(; approx_mode = "growth")
    sys = system(; approx_mode = approx_mode)
    return PB.SafetyProblem(sys, sys.X, sys.X, PB.Infinity())
end

end
