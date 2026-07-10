module DCDC

using StaticArrays
using MathematicalSystems
using Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem

Base.@kwdef struct Params{T}
    vs::T = 1.0
    xL::T = 3.0
    xC::T = 70.0
    r0::T = 1.0
    rL::T = 0.05
    rC::T = 0.005
end

function A1(p::Params = Params())
    return SMatrix{2, 2}(-p.rL / p.xL, 0.0, 0.0, -1.0 / (p.xC * (p.r0 + p.rC)))
end

function A2(p::Params = Params())
    return SMatrix{2, 2}(
        -(p.rL + p.r0 * p.rC / (p.r0 + p.rC)) / p.xL,
        5.0 * p.r0 / ((p.r0 + p.rC) * p.xC),
        -p.r0 / ((p.r0 + p.rC) * p.xL * 5.0),
        -1.0 / (p.xC * (p.r0 + p.rC)),
    )
end

function A2_abs(p::Params = Params())
    return SMatrix{2, 2}(
        -(p.rL + p.r0 * p.rC / (p.r0 + p.rC)) / p.xL,
        5.0 * p.r0 / ((p.r0 + p.rC) * p.xC),
        p.r0 / ((p.r0 + p.rC) * p.xL * 5.0),
        -1.0 / (p.xC * (p.r0 + p.rC)),
    )
end

function jacobian_bound(p::Params = Params())
    A1p = A1(p)
    A2p_abs = A2_abs(p)

    return u -> u[1] == 1 ? A1p : A2p_abs
end

function DF_sys(p::Params = Params())
    A1p = A1(p)
    A2p = A2(p)

    return u -> u[1] == 1 ? A1p : A2p
end

function dynamic(p::Params = Params())
    b = SVector(p.vs / p.xL, 0.0)
    A1p = A1(p)
    A2p = A2(p)

    return (x, u) -> u[1] == 1 ? A1p * x + b : A2p * x + b
end

function system(;
    params::Params = Params(),
    _X_ = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85)),
    _U_ = UT.HyperRectangle(SVector(1), SVector(2)),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        UT.get_dim(_X_),
        UT.get_dim(_U_),
        _X_,
        _U_,
    )
end

function problem(;
    params::Params = Params(),
    _I_ = UT.HyperRectangle(SVector(1.19, 5.59), SVector(1.21, 5.61)),
)
    sys = system(; params = params)

    return PR.SafetyProblem(sys, _I_, sys.X)
end

function system_plot!(;
    params::Params = Params(),
    xlims = (-0.5, 6.5),
    ylims = (-2.0, 2.0),
    show_state = true,
    show_params = false,
)
    return function (fig, x, u)
        iL = Float64(x[1])
        vC = Float64(x[2])
        mode = Int(round(u[1]))

        switch_closed = mode == 1

        plot!(fig, [0.0, 1.0], [0.0, 0.0]; lw = 3, color = :black, label = "")
        plot!(fig, [2.0, 3.0], [0.0, 0.0]; lw = 3, color = :black, label = "")
        plot!(fig, [4.0, 6.0], [0.0, 0.0]; lw = 3, color = :black, label = "")
        plot!(fig, [6.0, 6.0], [0.0, -1.0]; lw = 3, color = :black, label = "")
        plot!(fig, [0.0, 6.0], [-1.0, -1.0]; lw = 3, color = :black, label = "")
        plot!(fig, [0.0, 0.0], [-1.0, 0.0]; lw = 3, color = :black, label = "")

        scatter!(fig, [0.0], [-0.5]; markersize = 14, color = :white, label = "")
        annotate!(fig, -0.25, -0.5, text("Vs", 10))

        plot!(fig, [1.0, 2.0], [0.0, 0.0]; lw = 2, color = :black, label = "")
        annotate!(fig, 1.5, 0.35, text("L", 11))

        if switch_closed
            plot!(fig, [3.0, 4.0], [0.0, 0.0]; lw = 4, color = :green, label = "")
            annotate!(fig, 3.5, 0.35, text("S closed", 9))
        else
            plot!(fig, [3.0, 3.45], [0.0, 0.0]; lw = 3, color = :black, label = "")
            plot!(fig, [3.55, 4.0], [0.25, 0.0]; lw = 3, color = :red, label = "")
            annotate!(fig, 3.5, 0.55, text("S open", 9))
        end

        plot!(fig, [3.5, 3.5], [0.0, -1.0]; lw = 2, color = :black, label = "")
        annotate!(fig, 3.75, -0.5, text("D", 10))

        plot!(fig, [5.1, 5.1], [0.0, -1.0]; lw = 2, color = :black, label = "")
        plot!(fig, [4.9, 4.9], [0.0, -1.0]; lw = 2, color = :black, label = "")
        annotate!(fig, 5.25, -0.45, text("C", 11))

        plot!(fig, [6.0, 6.0], [0.0, -1.0]; lw = 2, color = :black, label = "")
        annotate!(fig, 6.25, -0.45, text("R", 11))

        if show_state
            annotate!(
                fig,
                0.2,
                1.45,
                text(
                    "mode = $mode\n" *
                    "iL = $(round(iL; digits = 3))\n" *
                    "vC = $(round(vC; digits = 3))",
                    10,
                ),
            )
        end

        if show_params
            annotate!(
                fig,
                4.2,
                1.45,
                text(
                    "Vs = $(params.vs)\n" *
                    "L = $(params.xL)\n" *
                    "C = $(params.xC)\n" *
                    "R = $(params.r0)",
                    9,
                ),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end
