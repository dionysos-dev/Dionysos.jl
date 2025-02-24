"""
The structure

    ControlSystem{N, T}

is the abstract type that defines a control system.
"""
abstract type ControlSystem{N, T} end

function get_f_eval(sys)
    return sys.f_eval
end

function RungeKutta4(F, x, u, tstep, nsub::Int)
    τ = tstep / nsub
    for _ in 1:nsub
        Fx1 = F(x, u)
        xrk = x + Fx1 * (τ / 2.0)
        Fx2 = F(xrk, u)
        xrk = x + Fx2 * (τ / 2.0)
        Fx3 = F(xrk, u)
        xrk = x + Fx3 * τ
        Fx4 = F(xrk, u)
        x = x + (Fx1 + Fx2 * 2.0 + Fx3 * 2.0 + Fx4) * (τ / 6.0)
    end
    return x
end

"""
    ControlSystemGrowth{N, T, F1 <: Function, F2 <: Function, F3 <: Function} <: ControlSystem{N, T}

is one implementation of the `ControlSystem` type for which we have a growth bound function.
"""
struct ControlSystemGrowth{N, T, F1 <: Function, F2 <: Function, F3 <: Function} <:
       ControlSystem{N, T}
    tstep::Float64
    sysnoise::SVector{N, T}
    measnoise::SVector{N, T}
    sys_map::F1
    growthbound_map::F2
    sys_inv_map::F3
end

function discretize_system_with_growth_bound(
    tstep,
    F_sys,
    L_growthbound,
    sysnoise::SVector{N, T},
    measnoise::SVector{N, T},
    nsys,
    ngrowthbound,
) where {N, T}
    sys_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N, T}
    end
    sys_inv_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N, T}
    end
    F_growthbound = let sysnoise = sysnoise
        (r, u) -> L_growthbound(u) * r + sysnoise
    end
    growthbound_map = let ngrowthbound = ngrowthbound
        (r::SVector{N, T}, u, tstep) ->
            RungeKutta4(F_growthbound, r, u, tstep, ngrowthbound)::SVector{N, T}
    end
    return ControlSystemGrowth(
        tstep,
        sysnoise,
        measnoise,
        sys_map,
        growthbound_map,
        sys_inv_map,
    )
end

"""
    ControlSystemLinearized{N, T, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function, } <: ControlSystem{N, T}

is one implementation of the `ControlSystem` type for which we have linearized the system map.
"""
struct ControlSystemLinearized{
    N,
    T,
    F1 <: Function,
    F2 <: Function,
    F3 <: Function,
    F4 <: Function,
} <: ControlSystem{N, T}
    tstep::Float64
    measnoise::SVector{N, T}
    sys_map::F1
    linsys_map::F2
    error_map::F3
    sys_inv_map::F4
end

function RungeKutta4Linearized(F, DF, x, dx, u, tstep, nsub::Int)
    τ = tstep / nsub
    for i in 1:nsub
        Fx1 = F(x, u)
        DFx1 = DF(x, u) * dx
        xrk = x + Fx1 * (τ / 2.0)
        dxrk = dx + DFx1 * (τ / 2.0)
        Fx2 = F(xrk, u)
        DFx2 = DF(xrk, u) * dxrk
        xrk = x + Fx2 * (τ / 2.0)
        dxrk = dx + DFx2 * (τ / 2.0)
        Fx3 = F(xrk, u)
        DFx3 = DF(xrk, u) * dxrk
        xrk = x + Fx3 * τ
        dxrk = dx + DFx3 * τ
        Fx4 = F(xrk, u)
        DFx4 = DF(xrk, u) * dxrk
        x += (Fx1 + Fx2 * 2.0 + Fx3 * 2.0 + Fx4) * (τ / 6.0)
        dx += (DFx1 + DFx2 * 2.0 + DFx3 * 2.0 + DFx4) * (τ / 6.0)
    end
    return (x, dx)
end

# Give an upper-bound on x(t) staisfying x'(t) ≦ a*x(t) + b*exp(2at)
function BoundSecondOrder(a, b, tstep)
    if a ≈ 0.0
        return b * tstep
    else
        ρ = exp(a * tstep)
        return (b / a) * ρ * (ρ - 1.0) / 2.0
    end
end

function discretize_system_with_linearization(
    tstep,
    F_sys,
    DF_sys,
    bound_DF,
    bound_DDF,
    measnoise::SVector{N, T},
    nsys,
) where {N, T}
    sys_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N, T}
    end
    linsys_map = let nsys = nsys
        (x::SVector{N, T}, dx::SMatrix{N, N, T}, u, tstep) -> RungeKutta4Linearized(
            F_sys,
            DF_sys,
            x,
            dx,
            u,
            tstep,
            nsys,
        )::Tuple{SVector{N, T}, SMatrix{N, N, T}}
    end
    sys_inv_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) ->
            RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N, T}
    end
    error_map =
        (r::T, u, tstep) -> BoundSecondOrder(bound_DF(u), bound_DDF(u), tstep) * (r^2)::T
    return ControlSystemLinearized(
        tstep,
        measnoise,
        sys_map,
        linsys_map,
        error_map,
        sys_inv_map,
    )
end

"""
    SimpleSystem{N, T, F <: Function, F2} <: ControlSystem{N, T}   

is one implementation of the `ControlSystem` type. 
"""
struct SimpleSystem{N, T, F <: Function, F2} <: ControlSystem{N, T}
    tstep::Float64
    measnoise::SVector{N, T}
    sys_map::F
    f::F2
end

function NewSimpleSystem(tstep, F_sys, measnoise::SVector{N, T}, nsys) where {N, T}
    sys_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N, T}
    end
    return SimpleSystem(tstep, measnoise, sys_map, F_sys)
end

"""
    EllipsoidalAffineApproximatedSystem{}   
    
is a system whose dynamics is a noisy constrained affine control discrete system whose cells are ellipsoids, with a bound on the Lipschitz constant.
"""
struct EllipsoidalAffineApproximatedSystem{}
    dynamics::Dict{UT.Ellipsoid, MS.NoisyConstrainedAffineControlDiscreteSystem}
    L::Dict{UT.Ellipsoid, Float64} # smoothness constant to bound error
end

function _getLipschitzConstants(J, xi, rules)
    L = zeros(Base.length(xi))
    for (i, g) in enumerate(eachrow(J))
        Hg_s = Symbolics.jacobian(g, xi) #gets symbolic hessian of the i-th component of f(x,u,w)
        Hg = Symbolics.substitute(Hg_s, rules)

        # f_aux = eval(build_function(Hg)[1]);
        # mat = Base.invokelatest(f_aux)
        mat = Symbolics.value.(Hg)
        #println(mat)
        if any(x -> isa(x, IntervalArithmetic.Interval), mat)
            mat = IntervalArithmetic.Interval.(mat)
            eigenbox = IntervalLinearAlgebra.eigenbox(mat)
            L[i] = abs(eigenbox).hi
        else
            L[i] = max(abs.(eigen(mat).values)...)
        end
    end
    return L
end

# MS.NoisyAffineControlDiscreteSystem
# is a system of the form x(k+1) = Ax(k)+Bu(k)+c+Ew(k) with x(k)∈X, u(k)∈U, w(k)∈W
function buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)
    n = Base.length(x)
    m = Base.length(u)
    p = Base.length(w)
    xi = vcat(x, u, w)
    x̄i = vcat(x̄, ū, w̄)
    Xi = X × U × W
    sub_rules_Xi = Dict(xi[i] => Xi[i] for i in 1:(n + m + p))

    Jx = Symbolics.jacobian(f, x)
    Ju = Symbolics.jacobian(f, u)
    Jw = Symbolics.jacobian(f, w)
    Jxi = hcat(Jx, Ju, Jw)

    L = _getLipschitzConstants(Jxi, xi, sub_rules_Xi)

    sub_rules_x̄i = Dict(xi[i] => x̄i[i] for i in 1:(n + m + p))
    function evalSym(x)
        # When x is a Vector{SymbolicUtils.BasicSymbolic{Real}}, 
        # one needs to substitute each element of the vector
        if eltype(x) <: Symbolics.SymbolicUtils.BasicSymbolic{Real}
            return [Symbolics.substitute(elem, sub_rules_x̄i) for elem in x]
        end
        return Float64.(Symbolics.value.(Symbolics.substitute(x, sub_rules_x̄i)))
    end

    A = evalSym(Jx)
    B = evalSym(Ju)
    E = evalSym(Jw)
    c = vec(evalSym(f) - A * x̄ - B * ū - E * w̄)
    return (MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W), L)
end

############################################
############################################
############################################
# Affine discrete time system approximaton with its Lypschitz constant and the region of validity
# of the Lipschitz constant

struct AffineApproximationDiscreteSystem #<: ControlSystem
    constrainedAffineSys::MS.NoisyConstrainedAffineControlDiscreteSystem
    L::Any
    f_eval::Any
    function AffineApproximationDiscreteSystem(sys, L)
        f_eval_fun(x, u, w) = sys.A * x + sys.B * u + sys.D * w + sys.c
        return new(sys, L, f_eval_fun)
    end
end
function AffineApproximationDiscreteSystem(A, B, c, E, X, U, W, L)
    contSys = MS.NoisyConstrainedAffineControlDiscreteSystem(A, B, c, E, X, U, W)
    return AffineApproximationDiscreteSystem(contSys, L)
end

############################################
############################################
############################################

struct SymbolicSystem
    fsymbolicT::Any
    fsymbolic::Any
    Ts::Any
    nx::Any
    nu::Any
    nw::Any
    x::Any
    u::Any
    w::Any
    ΔX::Any
    ΔU::Any
    ΔW::Any
    X::Any
    U::Any
    W::Any
    obstacles::Any
    f_eval::Any
    f_backward_eval::Any
    Uformat::Any
    Wformat::Any
end
