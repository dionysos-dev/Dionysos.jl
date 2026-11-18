# Accuracy of the classes of `ground_cells_common.jl` against the true range of the foot height, for
# the former rule of `RP.infeasible_cells` (`:margin`) and its four bounds (`RP.FOOT_BOUNDS`: the
# linear range, the linear range plus the worst-case second-order remainder, the per-cell second order,
# the corner values).
#
# The true extremes of `z` over a cell are searched, cell by cell, by a multi-start local search that is
# exact along each axis: `z`, as a function of one joint angle `qᵢ` (the others fixed), is a sinusoid
# `C + A cos t + B sin t` (a joint rotates the rest of the chain rigidly), whose three coefficients
# follow from three evaluations of the kinematics at `qᵢ, qᵢ ± π/2`, and whose extremes over an interval
# are at its ends or at `atan(B, A)`. Cyclic coordinate ascent from 81 starting points of the cell (the
# points of its 3 × 3 × 3 × 3 grid, corners included) gives the extremes. The search finds a point of the
# cell, so it can only underestimate the range: a class contradicted by it is certainly wrong, a class
# that survives it is very probably right (not a proof; the proof of the second-order ranges is the
# Taylor bound).
#
# For every cell classified `:below` (removed) or `:above` (kept), the class is "contradicted" if the
# extremes found reach the other side: `z_max ≥ 0` for a `:below` cell, `z_min < 0` for an `:above`
# one. Also printed: the classes of the true range (`:below` if `z_max < 0`, ...), which no range can
# beat, and the largest height reached by a `:below` cell / the lowest by an `:above` one. The results
# are saved with `Serialization` in `ground_cells_plots/extremes_<name>.jls` (read by
# `compare_ground_cells.jl`).
#
# Set `CELL_SIZE` and `STATE_BOUNDS = (q_lb, q_ub)` before including this file (the step box and
# 0.1 rad by default, as in `ground_cells_classify.jl`); `NAME` names the output file.
#
#   julia -t auto --project=problems/BipedRobot tests/check_ground_cells_extremes.jl

using Serialization

include(joinpath(@__DIR__, "ground_cells_common.jl"))

cell_size = @isdefined(CELL_SIZE) ? CELL_SIZE : SVector(0.1, 0.1, 0.1, 0.1)
q_lb, q_ub =
    @isdefined(STATE_BOUNDS) ? STATE_BOUNDS :
    (SVector(-0.5, 0.0, -0.8, -1.0), SVector(0.5, 1.0, 0.8, 0.0))
geometry = RP.default_geometry()

name = @isdefined(NAME) ? NAME : "c010_box"
extremes_file =
    @isdefined(EXTREMES_FILE) ? EXTREMES_FILE :
    joinpath(@__DIR__, "ground_cells_plots", "extremes_$name.jls")

# The maximum of `sign ⋅ z` over the cell `q ± r`, from the starts `q + r ⋅ s`, `s ∈ {-1, 0, 1}⁴`, by
# cyclic coordinate ascent with the exact maximization along each axis.
function extremum(q, r, sign, geometry; sweeps = 5)
    z(p) = RP.swing_foot_position(p, geometry)[2]
    best = -Inf
    for s in Iterators.product(ntuple(_ -> (-1.0, 0.0, 1.0), 4)...)
        p = q .+ r .* SVector(s)
        for _ in 1:sweeps
            for i in 1:4
                e = SVector(ntuple(a -> a == i ? π / 2 : 0.0, 4))
                up, down, here = z(p + e), z(p - e), z(p)
                B = (up - down) / 2
                A = here - (up + down) / 2
                a, b = sign * A, sign * B
                # `t` ranges over the interval that keeps `pᵢ` in the cell.
                tlo, thi = q[i] - r[i] - p[i], q[i] + r[i] - p[i]
                t = a * cos(tlo) + b * sin(tlo) >= a * cos(thi) + b * sin(thi) ? tlo : thi
                tstar = atan(b, a)
                tlo <= tstar <= thi && (t = tstar)
                p = setindex(p, p[i] + t, i)
            end
        end
        best = max(best, sign * z(p))
    end
    return sign * best
end

# The labels of the presentation and the ground bounds.
methods = [
    ("margin", :margin),
    ("first_order", :linear),
    ("second_order", :second_order_global),
    ("second_order_cell", :second_order_cell),
    ("corners", :corners),
]
classes = Dict(
    m => ground_cells_by_bound(q_lb, q_ub, cell_size, geometry; foot_bound = bound) for
    (m, bound) in methods
)
class0 = classes["margin"]
ncells = length(class0.class)
r = cell_size ./ 2
zmin = zeros(ncells)
zmax = zeros(ncells)
println(
    "True extremes of z over ",
    ncells,
    " cells (multi-start search, ",
    Threads.nthreads(),
    " threads)…",
)
elapsed = @elapsed Threads.@threads :static for idx in 1:ncells
    I = CartesianIndices(class0.class)[idx]
    pos = Tuple(I) .+ class0.kmin .- 1
    q = SVector{4}(MP.get_coord_by_pos(class0.q_grid, pos))
    zmin[idx] = extremum(q, r, -1, geometry)
    zmax[idx] = extremum(q, r, +1, geometry)
end
println("  done in ", round(elapsed; digits = 1), " s")

# The class of every cell under the true range, and the comparison of every method with it.
true_below = count(<(0), zmax)
true_above = count(>=(0), zmin)
println(
    "\ntrue range: below $true_below  above $true_above  across ",
    ncells - true_below - true_above,
)
results = Dict{String, Any}()
for (m, _) in methods
    code = vec(classes[m].class)
    below = findall(==(0x01), code)
    above = findall(==(0x03), code)
    res = (;
        below = length(below),
        above = length(above),
        across = ncells - length(below) - length(above),
        below_contradicted = count(i -> zmax[i] >= 0, below),
        above_contradicted = count(i -> zmin[i] < 0, above),
        highest_below = isempty(below) ? -Inf : maximum(zmax[below]),
        lowest_above = isempty(above) ? Inf : minimum(zmin[above]),
    )
    results[m] = res
    println(
        rpad(m, 19),
        ":below ",
        lpad(res.below, 7),
        " (contradicted ",
        lpad(res.below_contradicted, 4),
        ", highest z ",
        round(100 * res.highest_below; digits = 2),
        " cm);  :above ",
        lpad(res.above, 7),
        " (contradicted ",
        lpad(res.above_contradicted, 4),
        ", lowest z ",
        round(100 * res.lowest_above; digits = 3),
        " cm);  :across ",
        res.across,
    )
end
mkpath(dirname(extremes_file))
Serialization.serialize(
    extremes_file,
    (;
        results,
        true_range = (;
            below = true_below,
            above = true_above,
            across = ncells - true_below - true_above,
        ),
        ncells,
        cell_size = collect(cell_size),
        q_lb = collect(q_lb),
        q_ub = collect(q_ub),
        elapsed,
    ),
)
println("Saved to ", extremes_file)
