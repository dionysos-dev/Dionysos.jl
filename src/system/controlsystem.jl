abstract type ControlSystem{N, T} end

struct SimpleSystem{N, T, F <: Function, F2} <: ControlSystem{N, T}
    tstep::Float64
    measnoise::SVector{N, T}
    sys_map::F
    f::F2
end

function NewSimpleSystem(tstep, F_sys, measnoise::SVector{N, T}, nsys) where {N, T}
    sys_map = let nsys = nsys
        (x::SVector{N, T}, u, tstep) ->
            runge_kutta4(F_sys, x, u, tstep, nsys)::SVector{N, T}
    end
    return SimpleSystem(tstep, measnoise, sys_map, F_sys)
end

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
        if any(x -> isa(x, IA.Interval), mat)
            mat = IA.Interval.(mat)
            eigenbox = IL.eigenbox(mat)
            L[i] = abs(eigenbox).hi
        else
            L[i] = max(abs.(eigen(mat).values)...)
        end
    end
    return L
end

function buildAffineApproximation(f, x, u, w, x̄, ū, w̄, X, U, W)
    n = Base.length(x)
    m = Base.length(u)
    p = Base.length(w)
    xi = vcat(x, u, w)
    x̄i = vcat(x̄, ū, w̄)
    Xi = vcat(collect(X), collect(U), collect(W))
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
