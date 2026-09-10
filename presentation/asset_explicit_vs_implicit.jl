# Explicit versus implicit state mapping, measured on the path-planning benchmark.
#
# The mapping is the bijection between an integer state label and a grid cell.
#   explicit: a Dict pos -> id and a Vector id -> pos, materialised cell by cell
#   implicit: no table at all, the label IS the row-major index of the cell in a
#             rectangular index box, so id <-> pos is arithmetic
#
# Same grid, same problem, same answer. What changes is memory and lookup.
#
# Run: julia --project=test presentation/asset_explicit_vs_implicit.jl

ENV["GKSwstype"] = "100"

using StaticArrays, JuMP
import LazySets
using Symbolics, MathOptSymbolicAD
using Dionysos
const DI = Dionysos
const MP = DI.Mapping
const SY = DI.Symbolic

const WALLS =
    [([1.0, 0.0], [1.2, 9.0]), ([2.2, 0.0], [2.4, 5.0]), ([2.2, 6.0], [2.4, 10.0])]

function jacobian_bound(u)
    β = abs(u[1] / cos(atan(tan(u[2]) / 2)))
    return SMatrix{3, 3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, β, β, 0.0)
end

function run(; implicit::Bool)
    model = Model(Dionysos.Optimizer)
    x_low, x_upp = [0.0, 0.0, -pi - 0.4], [4.0, 10.0, pi + 0.4]
    @variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = [0.4, 0.4, 0.0][i])
    @variable(model, -1 <= u[1:2] <= 1)
    @expression(model, α, atan(tan(u[2]) / 2))
    @constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
    @constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
    @constraint(model, ∂(x[3]) == u[1] * tan(u[2]))
    @constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
    @constraint(model, final(x[2]) in MOI.Interval(0.3, 0.8))
    for (lo, hi) in WALLS
        @constraint(model, x[1:2] ∉ MOI.HyperRectangle(lo, hi))
    end

    set_attribute(model, "jacobian_bound", jacobian_bound)
    set_attribute(model, "time_step", 0.3)
    set_attribute(
        model,
        "state_grid",
        MP.GridFree(SVector(0.0, 0.0, 0.0), SVector(0.2, 0.2, 0.2)),
    )
    set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.3, 0.3)))
    set_attribute(model, "print_level", 0)

    if implicit
        set_attribute(model, "use_implicit_mapping", true)
        set_attribute(
            model,
            "mapping_region",
            LazySets.Hyperrectangle(; low = x_low, high = x_upp),
        )
    end

    total = @elapsed optimize!(model)
    abs_sys = get_attribute(model, "abstract_system")
    mapping = SY.get_state_mapping(abs_sys)

    return (;
        status = termination_status(model),
        total,
        n = MP.get_n_state(mapping),
        bytes = Base.summarysize(mapping),
        mapping = typeof(mapping).name.name,
    )
end

println("— explicit mapping —")
e = run(; implicit = false)
println("  ", e)

println("— implicit mapping —")
i = run(; implicit = true)
println("  ", i)

println()
println("mapping type      explicit : $(e.mapping)")
println("mapping type      implicit : $(i.mapping)")
println("labels            explicit : $(e.n)")
println("labels            implicit : $(i.n)")
println("mapping memory    explicit : $(round(e.bytes / 1024^2; digits = 2)) MB")
println(
    "mapping memory    implicit : $(round(i.bytes / 1024; digits = 2)) kB  " *
    "(÷$(round(e.bytes / i.bytes; digits = 0)))",
)
println("end-to-end        explicit : $(round(e.total; digits = 1)) s")
println("end-to-end        implicit : $(round(i.total; digits = 1)) s")
