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

# data-driven plot
function plot_transitions!(set, f_eval, c_eval, W; dims = [1, 2], N = 100)
    samples = UT.sample(set; N = N)
    nx = UT.get_dims(set)
    for x in samples
        unew = c_eval(x)
        wnew = UT.sample(W)
        xnew = f_eval(x, unew, wnew)
        plot!(UT.DrawArrow(SVector{nx}(x[dims]), SVector{nx}(xnew[dims])); color = :black)
    end
end

# data-driven plot
function plot_check_feasibility!(set1, set2, f_eval, c_eval, W; dims = [1, 2], N = 100)
    plot!(set1; dims = dims, color = :green)
    plot!(set2; dims = dims, color = :red)
    return plot_transitions!(set1, f_eval, c_eval, W; dims = dims, N = N)
end

function plot_controller_cost!(
    set,
    c_eval,
    cost_eval;
    N = 10,
    scale = 0.0001,
    dims = [1, 2],
    color = :white,
    linewidth = 1,
)
    samples = UT.sample(set; N = N)
    costs = []
    for x in samples
        unew = c_eval(x)
        push!(costs, cost_eval(x, unew))
    end
    vmin = minimum(costs)
    vmax = maximum(costs)
    colorMap = UT.Colormap([vmin, vmax], Colors.colormap("Blues"))
    P = (1 / scale) * Matrix{Float64}(I(2))
    plot!(set; color = color, linealpha = 1.0, linewidth = linewidth, linecolor = :black)
    for (i, x) in enumerate(samples)
        plot!(UT.Ellipsoid(P, x[dims]); color = UT.get_color(colorMap, costs[i]), lw = 0)
    end
    return plot!(colorMap)
end
