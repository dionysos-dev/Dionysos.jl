# Tests for Multibody2D.jl. Run (needs Symbolics in the active environment) with
#   julia --project=<env with Symbolics> tests/test_multibody2d.jl
using Test
using LinearAlgebra
using Symbolics
include(joinpath(@__DIR__, "..", "Multibody2D.jl"))
using .Multibody2D

# Evaluate a symbolic expression/array at numeric values given as a Dict of variable => value.
# `fold = Val(true)` makes `substitute` evaluate to plain numbers.
_ev1(x, vals) = Float64(Symbolics.value(Symbolics.substitute(x, vals; fold = Val(true))))
ev(x::AbstractArray, vals) =
    reshape(Float64[_ev1(x[i], vals) for i in eachindex(x)], size(x))
ev(x, vals) = _ev1(x, vals)
rand_vals(vars) = Dict(v => 2 * rand() - 1 for v in vars)

# Planar serial chain of n links (revolute joints q_i, link i of length l_{i-1} from its parent).
function chain(n)
    sv(name, i) = Symbolics.variable(name, i)
    q, qd = [sv(:q, i) for i in 1:n], [sv(:qd, i) for i in 1:n]
    l, lc, m, I = ([sv(nm, i) for i in 1:n] for nm in (:l, :lc, :m, :I))
    g = sv(:g, 1)
    sys = MultibodySystem(q, qd; gravity = g)
    parent = sys.world
    for i in 1:n
        f = add_frame!(sys, Symbol(:link, i), parent, i == 1 ? [0, 0] : [l[i - 1], 0], q[i])
        add_point!(sys, Symbol(:body, i), f, [lc[i], 0]; mass = m[i], inertia = I[i])
        parent = f
    end
    return sys, [l; lc; m; I; g]
end

# Compare the normal form against the second-order dynamics D q̈ = h + Bτ at a random state.
function check_normal_form(sys, B, params)
    n = length(sys.q)
    nf = normal_form(sys, B)
    act, un = nf.actuated, nf.unactuated
    na = length(act)
    D, h, V = mass_matrix(sys), bias_vector(sys), potential_energy(sys)

    vals = Dict{Num, Float64}()
    for v in [nf.state; nf.input]
        vals[v] = 2rand() - 1
    end
    for v in params
        vals[v] = 0.5 + rand()
    end
    qdot = ev(nf.velocities, vals)
    full = copy(vals)
    for i in 1:n
        full[sys.qd[i]] = qdot[i]
    end
    Dn, hn = ev(D, full), ev(h, full)
    τ, Bn = ev(nf.input, vals), ev(B, vals)
    qdd = Dn \ (hn + Bn * τ)
    xdot = ev(nf.xdot, vals)

    # Ḋ and V̇ along q̇ by central differences (independent of the symbolic derivatives).
    ε = 1e-6
    shifted(s) =
        merge(full, Dict(sys.q[i] => full[sys.q[i]] + s * ε * qdot[i] for i in 1:n))
    Ddot = (ev(D, shifted(1)) - ev(D, shifted(-1))) / 2ε
    Vdot = (ev(V, shifted(1)) - ev(V, shifted(-1))) / 2ε

    approx(a, b) = isapprox(a, b; atol = 1e-5, rtol = 1e-5)
    pu = ev(nf.state[(n + na + 1):end], vals)
    @test approx((Dn * qdot)[un], pu)                          # p_u = (D q̇)_u
    @test approx(xdot[1:n], qdot)                              # q̇
    @test approx(xdot[(n + 1):(n + na)], qdd[act])             # v̇_a = q̈_a
    @test approx(xdot[(n + na + 1):end], (Ddot * qdot + Dn * qdd)[un])  # ṗ_u
    @test abs(qdot' * hn + qdot' * Ddot * qdot / 2 + Vdot) < 1e-4        # h conserves energy
    # ṗ_u does not depend on τ; neither does anything but the v̇_a rows of the gain.
    nu = length(un)
    @test all(iszero, nf.gain[1:n, :]) && all(iszero, nf.gain[(n + na + 1):end, :])
    for e in nf.xdot[(n + na + 1):end]
        @test isempty(intersect(Symbolics.get_variables(e), nf.input))
    end
    return nf
end

@testset "Multibody2D" begin
    @testset "double pendulum (closed-form reference)" begin
        @variables q[1:2] qd[1:2] l1 l2 lc1 lc2 m1 m2 I1 I2 g
        sys = MultibodySystem(collect(q), collect(qd); gravity = g)
        link1 = add_frame!(sys, :link1, sys.world, [0, 0], q[1])
        link2 = add_frame!(sys, :link2, link1, [l1, 0], q[2])
        add_point!(sys, :body1, link1, [lc1, 0]; mass = m1, inertia = I1)
        add_point!(sys, :body2, link2, [lc2, 0]; mass = m2, inertia = I2)
        tip = add_point!(sys, :tip, link2, [l2, 0])

        D = mass_matrix(sys)
        D11 = I1 + I2 + m1 * lc1^2 + m2 * (l1^2 + lc2^2 + 2 * l1 * lc2 * cos(q[2]))
        D12 = I2 + m2 * (lc2^2 + l1 * lc2 * cos(q[2]))
        D22 = I2 + m2 * lc2^2
        # Symbolic equality after full reduction.
        @test iszero(trig_reduce(D[1, 1] - D11))
        @test iszero(trig_reduce(D[1, 2] - D12))
        @test iszero(trig_reduce(D[2, 1] - D12))
        @test iszero(trig_reduce(D[2, 2] - D22))
        # ... and the reduced form is the compact one (a single cos(q2), no cos(q1)).
        @test !occursin("q[1]", string(D[1, 2])) && !occursin("q[1]", string(D[2, 2]))

        T = kinetic_energy(sys)
        @test iszero(
            trig_reduce(T - (D11 * qd[1]^2 / 2 + D12 * qd[1] * qd[2] + D22 * qd[2]^2 / 2)),
        )

        V = potential_energy(sys)
        V_ref = g * (m1 * lc1 * sin(q[1]) + m2 * (l1 * sin(q[1]) + lc2 * sin(q[1] + q[2])))
        @test iszero(trig_reduce(V - V_ref))

        # Positions / velocities in several frames.
        @test iszero(
            trig_reduce(position(tip)[1] - (l1 * cos(q[1]) + l2 * cos(q[1] + q[2]))),
        )
        @test iszero(
            trig_reduce(position(tip)[2] - (l1 * sin(q[1]) + l2 * sin(q[1] + q[2]))),
        )
        @test isequal.(position(tip, link2), [l2, 0]) |> all
        @test iszero(trig_reduce(position(tip, link1)[2] - l2 * sin(q[2])))
        # Velocity of a point rigid in link2, seen from link2, is zero ...
        @test all(iszero, velocity(sys, tip; relative_to = link2))
        # ... seen from link1 it is l2 q̇2 in the tangential direction of link2.
        v = velocity(sys, tip; relative_to = link1)
        @test iszero(trig_reduce(v[1] + l2 * sin(q[2]) * qd[2]))
        @test iszero(trig_reduce(v[2] - l2 * cos(q[2]) * qd[2]))
        # World velocity expressed in link2 axes has norm² = |world velocity|².
        vw = velocity(sys, tip)
        v2 = velocity(sys, tip; expressed_in = link2)
        @test iszero(trig_reduce(sum(x -> x^2, v2) - sum(x -> x^2, vw)))
    end

    @testset "simplify = false agrees numerically" begin
        @variables q[1:3] qd[1:3] a b c m1 m2 m3 I1 I2 I3
        sys = MultibodySystem(collect(q), collect(qd); gravity = 9.81)
        f1 = add_frame!(sys, :f1, sys.world, [a, 0], q[1] + 0.3)
        f2 = add_frame!(sys, :f2, f1, [b, c], q[2])
        f3 = add_frame!(sys, :f3, f2, [c, -b], q[3] - q[1])   # angle depends on other coords
        add_point!(sys, :b1, f1, [0.1, 0.2]; mass = m1, inertia = I1)
        add_point!(sys, :b2, f2, [b, 0]; mass = m2, inertia = I2, com = [0.3, -0.1])
        add_point!(sys, :b3, f3, [0, 0]; mass = m3, inertia = I3)
        vars = [sys.q; sys.qd; a; b; c; m1; m2; m3; I1; I2; I3]
        for _ in 1:3
            vals = rand_vals(vars)
            for v in (m1, m2, m3, I1, I2, I3)   # physical (positive) masses and inertias
                vals[v] = 0.5 + rand()
            end
            @test ev(mass_matrix(sys), vals) ≈ ev(mass_matrix(sys; simplify = false), vals)
            @test ev(kinetic_energy(sys), vals) ≈
                  ev(kinetic_energy(sys; simplify = false), vals)
            @test ev(potential_energy(sys), vals) ≈
                  ev(potential_energy(sys; simplify = false), vals)
            D = ev(mass_matrix(sys), vals)
            @test D ≈ D'
            @test all(>(0), eigvals(Symmetric(D)))   # positive definite for positive masses/inertias
        end
    end

    @testset "frame translation depending on q (cart-pole)" begin
        @variables x θ xd θd M m l g
        sys = MultibodySystem([x, θ], [xd, θd]; gravity = g)
        cart = add_frame!(sys, :cart, sys.world, [x, 0], 0)
        pole = add_frame!(sys, :pole, cart, [0, 0], θ)
        add_point!(sys, :cart_mass, cart, [0, 0]; mass = M)
        add_point!(sys, :bob, pole, [0, -l]; mass = m)   # point mass, pole hangs along -y
        D = mass_matrix(sys)
        @test iszero(trig_reduce(D[1, 1] - (M + m)))
        @test iszero(trig_reduce(D[1, 2] - m * l * cos(θ)))
        @test iszero(trig_reduce(D[2, 2] - m * l^2))
        @test iszero(trig_reduce(potential_energy(sys) + m * g * l * cos(θ)))
    end

    @testset "velocity vs finite differences (sibling frames)" begin
        @variables q[1:3] qd[1:3]
        sys = MultibodySystem(collect(q), collect(qd))
        root = add_frame!(sys, :root, sys.world, [q[1], 0.5], 0.2)
        la = add_frame!(sys, :la, root, [0.3, 0], q[2])
        lb = add_frame!(sys, :lb, root, [-0.3, 0], q[3])
        p = add_point!(sys, :p, la, [0.7, 0.1])
        for frame in (sys.world, root, la, lb), expr in (frame, lb)
            v = velocity(sys, p; relative_to = frame, expressed_in = expr)
            qv = rand(3)
            qdv = rand(3)
            h = 1e-6
            vals = Dict(zip([sys.q; sys.qd], [qv; qdv]))
            pos(qq) = ev(position(p, frame), Dict(zip(sys.q, qq)))
            _, a = relative_transform(frame, expr)
            Rm = ev(rotation_matrix(a), Dict(zip(sys.q, qv)))
            fd = (pos(qv .+ h * qdv) - pos(qv .- h * qdv)) / 2h
            @test ev(v, vals) ≈ Rm * fd atol = 1e-6
        end
    end

    @testset "api errors and bookkeeping" begin
        @variables q qd
        sys = MultibodySystem([q], [qd])
        f = add_frame!(sys, :f, sys.world, [1, 0], q)
        @test sys[:f] === f && sys[:world] === sys.world
        @test_throws ArgumentError add_frame!(sys, :f, sys.world, [0, 0], 0)      # duplicate name
        @test_throws ArgumentError add_frame!(sys, :g, sys.world, [0, 0, 0], 0)   # bad translation
        @test_throws ArgumentError add_point!(sys, :p, f, [0, 0]; inertia = 1)    # inertia w/o mass
        other = MultibodySystem([q], [qd])
        @test_throws ArgumentError add_frame!(sys, :h, other.world, [0, 0], 0)    # foreign frame
        @test_throws KeyError sys[:nope]
        @test_throws DimensionMismatch MultibodySystem([q], [qd, qd])
        # No bodies -> zero energies.
        @test all(iszero, mass_matrix(sys))
        @test iszero(potential_energy(sys))
    end

    @testset "normal form" begin
        @testset "acrobot (closed form)" begin
            sys, params = chain(2)
            q = sys.q
            nf = normal_form(sys, selection_matrix(2, [2]))
            @test nf.actuated == [2] && nf.unactuated == [1]
            @test length(nf.state) == 4 &&
                  all(isequal.(nf.state[1:3], [q[1], q[2], sys.qd[2]]))
            # ṗ₁ = ∂L/∂q₁ = -∂V/∂q₁ (D does not depend on q₁), no τ.
            V = potential_energy(sys)
            @test iszero(trig_reduce(nf.xdot[4] + Symbolics.derivative(V, q[1])))
            # Gain on v̇₂ is D₁₁ / det(D); everything else is zero.
            D = mass_matrix(sys)
            vals = Dict(v => 0.5 + rand() for v in [params; q])
            detD = ev(D[1, 1] * D[2, 2] - D[1, 2]^2, vals)
            @test ev(nf.gain[3, 1], vals) ≈ ev(D[1, 1], vals) / detD
            @test all(iszero, nf.gain[[1, 2, 4], :])
            for _ in 1:3
                check_normal_form(sys, selection_matrix(2, [2]), params)
            end
        end

        @testset "one / two / three unactuated, and fully actuated" begin
            for (n, actuated) in
                ((4, [2, 3, 4]), (3, [3]), (4, [4]), (2, [1, 2]), (4, [1, 3]))
                sys, params = chain(n)
                B = selection_matrix(n, actuated)
                for _ in 1:2
                    nf = check_normal_form(sys, B, params)
                    @test nf.actuated == actuated
                end
            end
        end

        @testset "inputs, momenta and B validation" begin
            sys, _ = chain(3)
            @variables u1 u2 π3
            nf = normal_form(
                sys,
                selection_matrix(3, [1, 2]);
                input = [u1, u2],
                momenta = [π3],
            )
            @test all(isequal.(nf.input, [u1, u2])) && isequal(nf.state[end], π3)
            B = selection_matrix(3, [3, 1])               # column order = input order
            @test B == [0 1; 0 0; 1 0]
            @test normal_form(sys, B).actuated == [1, 3]
            @test_throws DimensionMismatch normal_form(sys, zeros(2, 1))
            @test_throws DimensionMismatch normal_form(sys, B; momenta = [π3, π3])
            @test_throws DimensionMismatch normal_form(sys, B; input = [u1])
            @test_throws ArgumentError selection_matrix(3, [1, 1])
            @test_throws ArgumentError selection_matrix(3, [4])
            # Lagrangian is T - V.
            @test iszero(
                trig_reduce(lagrangian(sys) - kinetic_energy(sys) + potential_energy(sys)),
            )
        end
    end
end
