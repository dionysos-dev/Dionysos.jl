
module Control

using ReachabilityAnalysis, Plots, StaticArrays, IntervalArithmetic


function test_1()
    ω = 1.2
    T = 2*pi / ω

    function duffing!(du, u, p, t) #@taylorize
        local α = -1.0
        local β = 1.0
        local δ = 0.3
        local γ = 0.37

        x, v = u
        f = γ * cos(ω * t)

        # write the nonlinear differential equations defining the model
        du[1] = u[2]
        du[2] = -α*x - δ*v - β*x^3 + f
    end

    # set of initial states
    X0 = Hyperrectangle(low=[0.9, -0.1], high=[1.1, 0.1])

    # formulate the initial-value problem
    prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim=2)

    # solve using a Taylor model set representation
    sol = solve(prob, tspan=(0.0, 20*T), alg=TMJets())
    #println(sol)
    #fig = plot(aspect_ratio = 1,legend = false)
    # plot the flowpipe in state-space
    fig = plot(sol, vars=(1, 2), xlab="x", ylab="v", lw=0.5, color=:red)
    display(fig)

end
function test_2()
    # x' = v
    # v' = -4x
    A = [0 1; -4. 0]
    X0 = Singleton([1.0, 0.0])
    prob = @ivp(X' = AX, X(0) ∈ X0)

    f(ΔT) = solve(prob, tspan=(0.0, 5.0), alg=GLGM06(δ=ΔT))

    plot(f(0.3), vars=(0, 1), lab="ΔT=0.3", color=:yellow)
    plot!(f(0.1), vars=(0, 1), lab="ΔT=0.1", color=:lightblue)
    plot!(f(0.05), vars=(0, 1), xlab="time", ylab="x(t)", lab="ΔT=0.05", color=:green)

    dom = 0:0.01:5.0
    fig = plot!(dom, cos.(2.0 * dom), lab="Analytic", color=:magenta)
    display(fig)
end
function test_3()
    function lotka_volterra!(du, u, p, t)
        local α, β, γ, δ = 1.5, 1.0, 3.0, 1.0
        x, y = u
        du[1] = x * (α - β*y)
        du[2] = -y * (γ - δ*x)
    end
    X₀ = Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2])
    prob = @ivp(x' = lotka_volterra!(x), dim: 2, x(0) ∈ X₀)

    sol = solve(prob, T=5.0)
    setrep(sol)

    sol = overapproximate(sol, Zonotope)
    setrep(sol)

    fig = plot(sol, vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:lightblue, lab="Flowpipe")
    plot!(X₀, color=:orange, lab="Xo")
    display(fig)
end
function test_4()
    function duffing!(dx, x, p, t) #@taylorize
        # write the nonlinear differential equations defining the model
        dx[1] = 1.0 + zero(x[1])
        dx[2] = 1.0 + zero(x[2])
    end

    # set of initial states
    X0 = Hyperrectangle(low=[0.0, 0.0], high=[1.0, 1.0])

    # formulate the initial-value problem
    prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim=2)
    # solve using a Taylor model set representation
    tsep = 5.0
    sol_forward = solve(prob, tspan=(0.0, tsep), alg=TMJets())
    #solz1 = overapproximate(sol_forward, Zonotope);
    println(sol_forward)
    # formulate the initial-value problem
    #prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim=2)
    #sol_backward = solve(prob, tspan=(0.0, -10.0), alg=TMJets())


    # plot the flowpipe in state-space
    fig = plot(sol_forward, vars=(1, 2), xlab="x", ylab="y", lw=0.5, color=:red)
    display(fig)
    #fig2 = plot(sol_backward, vars=(1, 2), xlab="x", ylab="y", lw=0.5, color=:red)
    #display(fig2)
end
#test_1()

function f1(A)
    n = length(A)
    i = 1
    j = 2
    while j <= n && n-j>=j-2
        if A[j] == A[1]
            if j == 2
                return true
            end
            for k=1:j-2
                if A[i+k] != A[j+k]
                    break
                end
                if k == j-2
                    return true
                end
            end
        end
        j+=1
    end
    return false
end

#A = [1,2,3,1,4,1,2,3,1,4]
#A = [1,2,1,3,5,1,2,1,3,5,8]

#println(f1(A))
using Base.Threads

function test_thread()
    max_threads = nthreads()
    println(max_threads)
    L = [1,2,3,4]
    Threads.@threads for (i,a) in enumerate(L)
           println(Threads.threadid())
    end
end
#test_4()

function test_box()
    B = IntervalBox(1.0..3.0, 2.0..4.0)
    r = f(B)
    c = SVector(1.0,1.0)
    l = c+B
    println()
    println(B)
    println(r)
    println(l)
end
end # module
