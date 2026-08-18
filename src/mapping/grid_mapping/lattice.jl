# Inter-sample safety of grid abstractions.
#
# A grid abstraction validates, for each (cell, input), only the reach set *at*
# the sampling time τ — never the trajectory during [0, τ]. When the admissible
# set has holes (obstacles carved out of `X`), a step whose displacement exceeds
# one cell per axis can land past a thin excluded region without ever producing
# an invalid target cell ("obstacle jump").
#
# No-jump guarantee: if (i) every cell intersecting the excluded region is
# removed from the domain (OUTER carving) and (ii) each step moves the state by
# at most one cell per axis, then the straight segment between two consecutive
# states stays inside `source_cell ∪ target_cell` — two cells certified free —
# so the continuous trajectory cannot cross the excluded region and no
# intermediate-time check is needed. Bounding τ by the *obstacle width* instead
# is not sound: it prevents stepping fully over the obstacle but still allows
# grazing a corner between two free endpoints.

"""
    intersample_safe_time_step(grid, velocity_bound)

Largest sampling time `τ` such that a system whose velocity satisfies
`|ẋᵢ| ≤ velocity_bound[i]` moves by at most one cell of `grid` per axis per
step: `τ = minᵢ hᵢ / velocity_bound[i]`.

Together with OUTER obstacle carving, respecting this bound guarantees that the
continuous trajectory between consecutive samples stays inside the union of the
two (certified admissible) cells, closing the inter-sample soundness gap of
grid abstractions over domains with holes. Axes with a zero velocity bound are
unconstrained.
"""
function intersample_safe_time_step(grid::Grid, velocity_bound::AbstractVector)
    h = get_h(grid)
    length(velocity_bound) == length(h) ||
        error("velocity_bound must have one entry per grid axis")
    return minimum(h ./ abs.(velocity_bound))
end

"""
    cells_on_segment(grid, a, b)

Positions of every grid cell whose closure the straight segment from `a` to
`b` passes through, in traversal order. Boundary crossings are enumerated
exactly (one interval per pair of consecutive crossing times, quantized at the
interval midpoint), so each segment point lies in the closure of a returned
cell.

This is the sound validation for *multi-cell* steps of integrator-type
dynamics: their inter-sample trajectory is exactly this segment, so a
transition whose swept cells are all admissible cannot cross an excluded
region between samples — lifting the one-cell-per-step restriction of
[`intersample_safe_time_step`](@ref).
"""
function cells_on_segment(grid::Grid{N}, a, b) where {N}
    d = b .- a
    orig = get_origin(grid)
    h = get_h(grid)

    ts = [0.0, 1.0]
    for i in 1:N
        d[i] == 0 && continue
        lo, hi = minmax(a[i], b[i])
        # Cell faces of axis `i` lie at orig + (k + 1/2)·h.
        kmin = ceil(Int, (lo - orig[i]) / h[i] - 0.5)
        kmax = floor(Int, (hi - orig[i]) / h[i] + 0.5)
        for k in kmin:kmax
            face = orig[i] + (k + 0.5) * h[i]
            t = (face - a[i]) / d[i]
            0.0 < t < 1.0 && push!(ts, t)
        end
    end
    sort!(ts)

    cells = NTuple{N, Int}[]
    for j in 1:(length(ts) - 1)
        ts[j + 1] - ts[j] > 1e-12 || continue
        tm = (ts[j] + ts[j + 1]) / 2
        pos = get_pos_by_coord(grid, a .+ tm .* d)
        (isempty(cells) || cells[end] != pos) && push!(cells, pos)
    end
    return cells
end

"""
    is_lattice_exact(state_grid, input_grid, tstep; rtol = 1e-9)

For integrator dynamics `ẋ = u`, check that every input-grid point `u`
translates the state grid exactly onto itself in one step: `tstep .* u` must be
a per-axis integer multiple of the state-grid step (checked on the input grid's
step and origin, which spans all its points).

When this holds, the exact reach set of a cell is again a single cell, so the
abstraction is deterministic and *exact* (a bisimulation) — in particular the
`CENTER_SIMULATION` approximation mode, unsound in general, becomes exact for
such systems.
"""
function is_lattice_exact(state_grid::Grid, input_grid::Grid, tstep::Real; rtol = 1e-9)
    hx = get_h(state_grid)
    du = get_h(input_grid)
    u0 = get_origin(input_grid)
    length(du) == length(hx) || return false
    on_lattice(v, h) = begin
        q = v / h
        return abs(q - round(q)) <= rtol * max(1.0, abs(q))
    end
    return all(
        on_lattice(tstep * du[i], hx[i]) && on_lattice(tstep * u0[i], hx[i]) for
        i in eachindex(hx)
    )
end
