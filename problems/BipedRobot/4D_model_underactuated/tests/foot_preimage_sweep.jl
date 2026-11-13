# The pre-image of a point of the swing-foot plane, and the cells of each class that it crosses.
#
# Not a script: it only defines functions, shared by `plot_ground_preimage_counts.jl` (which explains
# the maths) and `ground_cells_linearized.jl`. Include it after defining `MP` (`Dionysos.Mapping`) and
# `RP` (the robot problem module), with `StaticArrays` loaded. The cells are given as an array indexed
# by position whose entries are class codes `1, 2, …, K` (see `cell_class_array`).

# Swing-leg inverse kinematics from `c = cos q₄` (known to be in [-1, 1]) and `d = foot - hip`, for the
# knee branch `sgn` (the sign of `q₄`): the absolute angle `φ₃ = q₁ + q₂ + q₃` of the swing femur and
# `q₄`. The formulas of `RP._swing_leg_ik`, which only does the branch with the knee bent backwards.
@inline function swing_leg_ik(dx, dz, c, sgn, L1, L2)
    q4 = sgn * acos(c)
    A = L1 + L2 * c
    B = sgn * L2 * sqrt(1 - c * c)
    return atan(dx, -dz) - atan(B, A), q4
end

# The cells as an array indexed by position, with the class code `codes[k]` (a `UInt8`, from 1) of
# the cell `positions[k]`. The cells of a box form a box of positions `kmin:kmax`, so every index of
# the array is a cell.
function cell_class_array(positions, codes)
    kmin = ntuple(a -> minimum(pos[a] for pos in positions), 4)
    kmax = ntuple(a -> maximum(pos[a] for pos in positions), 4)
    class = zeros(UInt8, (kmax .- kmin .+ 1)...)
    for (pos, code) in zip(positions, codes)
        class[(pos .- kmin .+ 1)...] = code
    end
    return class, kmin, kmax
end

# Everything the sweep needs, for the domain `kmin:kmax` of cells: the hip at the `(q₁, q₂)` samples
# (`step`, midpoint rule, so that no sample lies on a face between two cells), the cell index of the
# samples, and the range of `cos q₄` of each knee branch (`q₄ = ±acos(c)` lies in the domain).
function build_sweep(geometry, class, kmin, kmax, h, step)
    L1, L2 = geometry.femur, geometry.tibia
    lo = SVector{4}(ntuple(a -> (kmin[a] - 0.5) * h[a], 4))
    hi = SVector{4}(ntuple(a -> (kmax[a] + 0.5) * h[a], 4))
    n1 = ceil(Int, (hi[1] - lo[1]) / step)
    n2 = ceil(Int, (hi[2] - lo[2]) / step)
    q1 = [lo[1] + (i - 0.5) * (hi[1] - lo[1]) / n1 for i in 1:n1]
    q2 = [lo[2] + (j - 0.5) * (hi[2] - lo[2]) / n2 for j in 1:n2]
    hip = [
        RP.cartesian_coordinates(SVector(a, b, 0.0, 0.0), geometry)[3] for a in q1, b in q2
    ]
    dims = size(class)
    k1 = [round(Int, a / h[1]) for a in q1]
    k2 = [round(Int, b / h[2]) for b in q2]
    # acos(c) in [a, b] <=> c in [cos b, cos a]; |q₄| ≤ π, and q₄ = sgn * acos(c).
    branches = Tuple{Float64, Float64, Float64}[]
    for sgn in (1.0, -1.0)
        a, b = sgn > 0 ? (max(lo[4], 0.0), hi[4]) : (max(-hi[4], 0.0), -lo[4])
        a <= b && push!(branches, (sgn, cos(b), cos(a)))
    end
    return (;
        L1,
        L2,
        h,
        kmin,
        kmax,
        strides = (1, dims[1], dims[1] * dims[2], dims[1] * dims[2] * dims[3]),
        class,
        hx = [p[1] for p in hip],
        hz = [p[2] for p in hip],
        θ = [a + b for a in q1, b in q2],
        # Linear index (offset) of the cell of each (q₁, q₂) sample, without the (q₃, q₄) part.
        offset12 = [
            (k1[i] - kmin[1]) + (k2[j] - kmin[2]) * dims[1] for i in 1:n1, j in 1:n2
        ],
        branches,
    )
end

# The number of cells of each class crossed by the pre-image of the point `(px, pz)`, written into
# `counts` (one entry per class code). `seen` and `touched` are scratch space (one cell flag, and the
# flags to reset), left clean.
function preimage_counts!(
    counts::Vector{Int},
    seen::Vector{Bool},
    touched::Vector{Int},
    sw,
    px,
    pz,
)
    fill!(counts, 0)
    L1, L2 = sw.L1, sw.L2
    c0 = L1^2 + L2^2
    k = 2 * L1 * L2
    h3, h4 = sw.h[3], sw.h[4]
    s3, s4 = sw.strides[3], sw.strides[4]
    for (sgn, c_lo, c_hi) in sw.branches
        for j in axes(sw.hx, 2), i in axes(sw.hx, 1)
            dx = px - sw.hx[i, j]
            dz = pz - sw.hz[i, j]
            c = (dx * dx + dz * dz - c0) / k
            # The swing foot is reachable from this hip, with a knee angle in the domain.
            c_lo <= c <= c_hi || continue
            φ3, q4 = swing_leg_ik(dx, dz, c, sgn, L1, L2)
            k4 = round(Int, q4 / h4)
            sw.kmin[4] <= k4 <= sw.kmax[4] || continue
            q3 = mod(φ3 - sw.θ[i, j] + π, 2π) - π
            k3 = round(Int, q3 / h3)
            sw.kmin[3] <= k3 <= sw.kmax[3] || continue
            idx = 1 + sw.offset12[i, j] + (k3 - sw.kmin[3]) * s3 + (k4 - sw.kmin[4]) * s4
            seen[idx] && continue
            seen[idx] = true
            push!(touched, idx)
            counts[sw.class[idx]] += 1
        end
    end
    for idx in touched
        seen[idx] = false
    end
    empty!(touched)
    return counts
end

# A point of the pre-image of `p` above `(q₁, q₂)` on the knee branch `sgn`, or `nothing` if the foot
# cannot reach `p` from that hip. Uses the same inverse kinematics as the sweep.
function preimage_point(geometry, q1, q2, p, sgn)
    L1, L2 = geometry.femur, geometry.tibia
    hip = RP.cartesian_coordinates(SVector(q1, q2, 0.0, 0.0), geometry)[3]
    d = p - hip
    c = (d[1]^2 + d[2]^2 - L1^2 - L2^2) / (2 * L1 * L2)
    abs(c) <= 1 || return nothing
    φ3, q4 = swing_leg_ik(d[1], d[2], c, sgn, L1, L2)
    return SVector(q1, q2, φ3 - (q1 + q2), q4)
end

# The window of the plane that holds the image of the domain: the bounding box of the swing foot at
# the corners and the centre of every cell, plus a margin.
function foot_window(positions, q_grid, geometry, margin)
    h = MP.get_h(q_grid)
    xmin = zmin = Inf
    xmax = zmax = -Inf
    for pos in positions
        c = MP.get_coord_by_pos(q_grid, pos)
        points = [c]
        append!(
            points,
            [
                c .+ h .* SVector(s) for
                s in Iterators.product(ntuple(_ -> (-0.5, 0.5), 4)...)
            ],
        )
        for q in points
            foot = RP.swing_foot_position(SVector{4}(q), geometry)
            xmin, xmax = min(xmin, foot[1]), max(xmax, foot[1])
            zmin, zmax = min(zmin, foot[2]), max(zmax, foot[2])
        end
    end
    return xmin - margin, xmax + margin, zmin - margin, zmax + margin
end

# The counts of the `nclasses` classes at every pixel `(xs[i], zs[j])`, as an array `[i, j, class]`,
# in parallel over the rows.
function count_plane(sw, xs, zs, nclasses)
    nthreads = Threads.maxthreadid()
    seen = [fill(false, length(sw.class)) for _ in 1:nthreads]
    touched = [Int[] for _ in 1:nthreads]
    scratch = [zeros(Int, nclasses) for _ in 1:nthreads]
    counts = zeros(Int, length(xs), length(zs), nclasses)
    Threads.@threads :static for j in eachindex(zs)
        t = Threads.threadid()
        for i in eachindex(xs)
            preimage_counts!(scratch[t], seen[t], touched[t], sw, xs[i], zs[j])
            counts[i, j, :] .= scratch[t]
        end
    end
    return counts
end
