abstract type Controller end

function get_c_eval(cont)
    return cont.c_eval
end

# Constant state-dependent controller of the form: κ(x) = c 
struct ConstantController{T<:Real, VT<:AbstractVector{T}} <: Controller
    c::VT
    c_eval

    function ConstantController(c::VT) where {T<:Real, VT<:AbstractVector{T}}
        c_eval_fun(x) = c
        return new{T,VT}(c, c_eval_fun)
    end
end

# Affine state-dependent controller of the form: κ(x) = K*(x-c)+ℓ
struct AffineController{T<:Real, MT<:AbstractMatrix{T}, VT1<:AbstractVector{T}, VT2<:AbstractVector{T}} <: Controller
    K::MT
    c::VT1
    ℓ::VT2
    c_eval

    function AffineController(K::MT, c::VT1, ℓ::VT2) where {T<:Real, MT<:AbstractMatrix{T}, VT1<:AbstractVector{T}, VT2<:AbstractVector{T}}
        c_eval_fun(x) = K * (x - c) + ℓ
        return new{T,MT,VT1,VT2}(K, c, ℓ, c_eval_fun)
    end
end

# data-driven check
function check_feasibility(set1, set2, f_eval, c_eval, nw, Uset; N=500)
    samples = UT.sample(set1; N=N)
    wnew = zeros(nw)
    for x in samples
        unew = c_eval(x)
        if !(unew ∈ Uset)
            println("Not feasible input")
            return false
        end
        xnew = f_eval(x, unew, wnew)
        if !(xnew ∈ set2)
            println("Not in the target ellipsoid")
            return false
        end
    end
    return true
end

# data-driven plot
function plot_transitions!(set, f_eval, c_eval, nw; dims=[1,2], N=100)
    samples = UT.sample(set; N=N)
    nx = UT.get_dims(set)
    wnew = zeros(nw)
    for x in samples
        unew = c_eval(x)
        xnew = f_eval(x, unew, wnew)
        plot!(UT.DrawArrow(SVector{nx}(x[dims]), SVector{nx}(xnew[dims])), color = :black)
    end
end

# data-driven plot
function plot_check_feasibility!(set1, set2, f_eval, c_eval, nw; dims=[1,2], N=100)
    fig = plot(aspect_ratio=:equal)
    plot!(set1, dims=dims, color = :green)
    plot!(set2, dims=dims, color = :red)
    plot_transitions!(set1, f_eval, c_eval, nw; dims=dims, N=N)
    display(fig)
end

function plot_controller_cost!(set, c_eval, cost_eval; N=10, scale=0.0001, dims=[1,2], color=:white)
    samples = UT.sample(set; N=N)
    costs = []
    for x in samples
        unew = c_eval(x)  
        push!(costs, cost_eval(x, unew))
    end
    vmin = minimum(costs)
    vmax = maximum(costs)
    colorMap = UT.Colormap([vmin,vmax], Colors.colormap("Blues"))
    P = (1/scale)*Matrix{Float64}(I(2))
    plot!(set, color=color)
    for (i,x) in enumerate(samples)
        plot!(UT.Ellipsoid(P, x[dims]), color=UT.get_color(colorMap, costs[i]), lw=0)
    end
    plot!(colorMap)
end



