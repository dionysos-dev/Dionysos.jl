abstract type Controller end

function get_c_eval(cont)
    return cont.c_eval
end

# Constant state-dependent controller of the form: κ(x) = c 
"""
    ConstantController{T, VT}

encodes a constant state-dependent controller of the κ(x) = c.
"""
struct ConstantController{T <: Real, VT <: AbstractVector{T}} <: Controller
    c::VT
    c_eval::Any

    function ConstantController(c::VT) where {T <: Real, VT <: AbstractVector{T}}
        c_eval_fun(x) = c
        return new{T, VT}(c, c_eval_fun)
    end
end

# Affine state-dependent controller of the form: κ(x) = K*(x-c)+ℓ
"""
    AffineController{T, MT, VT1, VT2}

encodes an affine state-dependent controller of the κ(x) = K*(x-c)+ℓ.
"""
struct AffineController{
    T <: Real,
    MT <: AbstractMatrix{T},
    VT1 <: AbstractVector{T},
    VT2 <: AbstractVector{T},
} <: Controller
    K::MT
    c::VT1
    ℓ::VT2
    c_eval::Any

    function AffineController(
        K::MT,
        c::VT1,
        ℓ::VT2,
    ) where {
        T <: Real,
        MT <: AbstractMatrix{T},
        VT1 <: AbstractVector{T},
        VT2 <: AbstractVector{T},
    }
        c_eval_fun(x) = K * (x - c) + ℓ
        return new{T, MT, VT1, VT2}(K, c, ℓ, c_eval_fun)
    end
end

struct ControllerAnalysis{C, F, W, S1, S2, CE}
    c_eval::C
    f_eval::F
    Wset::W
    domain_set::S1
    target_set::Union{Nothing, S2}
    cost_eval::Union{Nothing, CE}
    dims::Vector{Int}
    N::Int
end

function ControllerAnalysis(
    c_eval,
    f_eval,
    Wset,
    domain_set;
    target_set = nothing,
    cost_eval = nothing,
    dims = [1, 2],
    N = 100,
)
    return ControllerAnalysis(
        c_eval,
        f_eval,
        Wset,
        domain_set,
        target_set,
        cost_eval,
        dims,
        N,
    )
end

@recipe function f(
    analysis::ControllerAnalysis;
    arrowsB = true,
    cost = true,
    sphere_radius = 0.01,
)
    dims = analysis.dims
    samples = UT.sample(analysis.domain_set; N = analysis.N)

    c_eval = analysis.c_eval
    nx = UT.get_dims(analysis.domain_set)
    # Show domain and target sets
    @series begin
        color := :green
        return analysis.domain_set
    end

    if analysis.target_set !== nothing
        @series begin
            color := :red
            return analysis.target_set
        end
    end

    # Compute cost (if enabled)
    if cost && analysis.cost_eval !== nothing
        costs = Float64[]
        for x in samples
            u = c_eval(x)
            push!(costs, analysis.cost_eval(x, u))
        end
        vmin, vmax = extrema(costs)
        mycolorMap = UT.Colormap([vmin, vmax], Colors.colormap("Blues"))
        P = (1.0 / sphere_radius) * Matrix{Float64}(I(2))  # small ellipsoid

        for (i, x) in enumerate(samples)
            color = UT.get_color(mycolorMap, costs[i])
            ell = UT.Ellipsoid(P, x[dims])
            @series begin
                color := color
                lw := 0
                return ell
            end
        end
        @series begin
            mycolorMap
        end
    end

    if arrowsB
        for x in samples
            u = c_eval(x)
            w = UT.sample(analysis.Wset)
            x2 = analysis.f_eval(x, u, w)

            @series begin
                color := :black
                return UT.DrawArrow(SVector{nx}(x[dims]), SVector{nx}(x2[dims]))
            end
        end
    end
end

# data-driven check
function check_feasibility(
    ell1,
    ell2,
    f_eval,
    c_eval,
    Uset,
    Wset;
    N = 500,
    input_check = true,
    noise_check = true,
)
    samples = UT.sample(ell1; N = N)
    nw = UT.get_dims(Wset)
    for x in samples
        unew = c_eval(x)
        if input_check && !(unew ∈ Uset)
            println("Not feasible input")
            return false
        end
        noise_check ? wnew = UT.sample(Wset) : wnew = zeros(nw)
        xnew = f_eval(x, unew, wnew)
        if !(xnew ∈ ell2)
            println("Not in the target ellipsoid")
            return false
        end
    end
    return true
end
