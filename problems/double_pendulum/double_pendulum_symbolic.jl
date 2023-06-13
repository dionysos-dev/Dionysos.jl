using RigidBodyDynamics, MeshCat, MeshCatMechanisms, SymPy #Symbolics
using LinearAlgebra, StaticArrays, Random

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const PB = DI.Problem
const ST = DI.System
const CO = DI.Control
const SY = DI.Symbolic

Random.seed!(123)

include("double_pendulum_urdf.jl")

# q1: is the shoulder joint angle
# q2: is the elbow (relative) joint angle
# v1 and v2 are the respective time derivative
# we use x = (q1,q2,v1,v2)
# The zero configuration is with both links pointed directly down.
# The moments of inertia are taken about the pivots.
# http://underactuated.csail.mit.edu/acrobot.html#section1
function f(
    x,
    u = [0.0];
    g = 9.81,
    I1 = 0.333,
    lc1 = 0.5,
    m1 = 1.0,
    l1 = 1.0,
    I2 = 1.33,
    lc2 = 1.0,
    m2 = 1.0,
)
    q1 = x[1]
    q2 = x[2]
    v1 = x[3]
    v2 = x[4]
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    s12 = sin(q1 + q2)

    M11 = I1 + I2 + m2 * l1^2 + 2 * m2 * l1 * lc2 * c2
    M12 = I2 + m2 * l1 * lc2 * c2
    M22 = I2
    M = @SMatrix [M11 M12; M12 M22]

    C11 = -2 * m2 * l1 * lc2 * s2 * v2
    C12 = -m2 * l1 * lc2 * s2 * v2
    C21 = m2 * l1 * lc2 * s2 * v1
    C22 = 0
    C = @SMatrix [C11 C12; C21 C22]

    G = @SMatrix [m1 * g * lc1 * s1 + m2 * g * (l1 * s1 + lc2 * s12); m2 * g * lc2 * s12]
    v = @SMatrix [v1; v2]
    τ = @SMatrix [u[1]; 0.0] # for pendubot @SMatrix [0.0;u[1]] for acrobot
    r = inv(M) * (-C * v - G + τ)
    F = [x[3], x[4], r[1], r[2]] # IntervalBox(x[3],x[4],r[1],r[2])
    return F
end

function NewControlSystemFunction(f)
    function sys_map(x, u, tstep)
        nsub = 4
        return ST.RungeKutta4(f, x, u, tstep, nsub)
    end
end

function NewControlSystemUrdf(urdf)
    mechanism = parse_urdf(urdf)
    function sys_map(x, u, tstep)
        q = [x[1], x[2]]
        v = [x[3], x[4]]
        function control!(torques::AbstractVector, t, state::MechanismState)
            torques[1] = u[1]
            return torques[2] = 0.0
        end
        state = MechanismState(mechanism, q, v)
        t, q, v = simulate(state, tstep, control!; Δt = tstep)
        return SVector(q[2]..., v[2]...)
    end
    return sys_map
end

function check_one_config(urdf_file, f)
    sys_map1 = NewControlSystemUrdf(urdf_file)
    sys_map2 = NewControlSystemFunction(f)

    x = [0.5; 0.4; 2.0; -1.0]
    u = [0.4]
    tstep = 0.1

    next1 = sys_map1(x, u, tstep)
    next2 = sys_map2(x, u, tstep)
    err = norm(next1 - next2)
    println(err)
    println(next1)
    return println(next2)
end

# note: when tstep decreases, the number of errors decreases.
# when we have errors, they come more from the time derivative of the angles
function check_dynamical_equations(urdf_file, f; K = 10, verbose = true)
    count = 0
    check = true
    tol = 0.2
    for k in 1:K
        x = rand(4) * 2π
        u = rand(1) * 0.4
        tstep = rand() * 0.06
        sys_map1 = NewControlSystemUrdf(urdf_file)
        sys_map2 = NewControlSystemFunction(f)
        next1 = sys_map1(x, u, tstep)
        next2 = sys_map2(x, u, tstep)
        err = norm(next1 - next2)
        if norm(err) > tol
            count = count + 1
            if verbose
                println(err)
                println(next1)
                println(next2)
            end
            check = false
        end
    end
    if verbose
        println("error: ", count, "/", K)
    end
    return check
end

urdf_file = "C:\\Users\\jcalbert\\Documents\\GitHub\\Dionysos.jl\\problems\\double_pendulum\\double_pendulum.urdf"
DoublePendulum.create_urdf_file(urdf_file)
# DoublePendulum.simulate_double_pendulum(urdf_file)

check_one_config(urdf_file, f)
println(check_dynamical_equations(urdf_file, f; K = 1000))
