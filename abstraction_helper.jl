function input_from_index(u_idx::Int, abstract_system)
    Udom = abstract_system.Udom

    # `Udom` is a `CustomList{2,Float64}` whose printed form shows
    # it wraps a `Vector{SVector{2,Float64}}`. We access that inner
    # vector via `getfield(Udom, 1)` to avoid depending on field names.
    u_list = getfield(Udom, 1)::Vector{SVector{2,Float64}}

    @assert 1 ≤ u_idx ≤ length(u_list) "u_idx=$(u_idx) out of range 1:$(length(u_list))"

    return u_list[u_idx]
end

function _grid_params(grid)
    # 1) Try common field/property names across versions
    candidates_x0 = Any[]
    candidates_hx = Any[]

    # property-based access (works if fields are public)
    for nm in (:x0, :origin, :lb, :lower, :min)
        if hasproperty(grid, nm)
            push!(candidates_x0, getproperty(grid, nm))
        end
    end
    for nm in (:h, :hx, :step, :Δ, :delta)
        if hasproperty(grid, nm)
            push!(candidates_hx, getproperty(grid, nm))
        end
    end

    # 2) If not found, scan actual struct fields and pick sensible SVectors/Tuples
    # (GridFree in some Dionysos versions stores these under different names)
    for fn in fieldnames(typeof(grid))
        v = try
            getfield(grid, fn)
        catch
            continue
        end
        # accept StaticArrays SVectors, Tuples, or small vectors
        if v isa StaticArrays.SVector || v isa NTuple || v isa Tuple || v isa AbstractVector
            # we only care about length-3 objects
            lv = try length(v) catch; 0 end
            if lv == 3
                # classify by magnitude: hx should be strictly positive and typically small
                allpos = try all(t -> t > 0, v) catch; false end
                if allpos
                    push!(candidates_hx, v)
                else
                    push!(candidates_x0, v)
                end
            end
        end
    end

    # 3) Final fallback: use the script globals if they exist
    if isempty(candidates_x0)
        if @isdefined(x0)
            push!(candidates_x0, x0)
        end
    end
    if isempty(candidates_hx)
        if @isdefined(hx)
            push!(candidates_hx, hx)
        end
    end

    x0v = isempty(candidates_x0) ? nothing : candidates_x0[1]
    hxv = isempty(candidates_hx) ? nothing : candidates_hx[1]

    if x0v === nothing || hxv === nothing
        # Emit a helpful, actionable error
        error(
            "Could not extract grid parameters (x0,hx) from grid type=$(typeof(grid)). " *
            "Try in the REPL: fieldnames(typeof(abstract_system.Xdom.grid)) and also: " *
            "map(f->(f, getfield(abstract_system.Xdom.grid,f)), fieldnames(typeof(abstract_system.Xdom.grid)))"
        )
    end

    # Normalize to SVector{3,Float64}
    x0s = x0v isa StaticArrays.SVector ? x0v : SVector{3,Float64}(Float64.(collect(x0v))...)
    hxs = hxv isa StaticArrays.SVector ? hxv : SVector{3,Float64}(Float64.(collect(hxv))...)

    return x0s, hxs
end

"""Return the integer grid position (i,j,k) for abstract state id s."""
function _pos_of_state(abs_sys, s::Int)
    if hasproperty(abs_sys, :XMapping)
        mapping = abs_sys.XMapping
        if hasproperty(mapping, :id2pos)
            return mapping.id2pos[s]
        end
    end

    # Backward compatibility with older Dionysos versions.
    for fname in (:xint2pos, :id2pos, :state_to_pos)
        if hasproperty(abs_sys, fname)
            mapping = getproperty(abs_sys, fname)
            if mapping isa AbstractVector || mapping isa AbstractDict
                return mapping[s]
            end
        end
    end

    error("Could not recover state-position mapping from abstract_system; fields=$(fieldnames(typeof(abs_sys)))")
end

"""Compute the 3D center of abstract state s from (xint2pos, Xdom.grid)."""
function _center3d_from_grid(abs_sys, s::Int)
    pos = _pos_of_state(abs_sys, s)
    grid = hasproperty(abs_sys, :XMapping) ? abs_sys.XMapping.grid : getproperty(getproperty(abs_sys, :Xdom), :grid)
    x0g, hxg = _grid_params(grid)

    # `pos` is typically a Tuple{Int,Int,Int}.
    # For GridFree, the representative point is at: x = x0 + (pos - 0.5) .* hx
    # This matches the earlier behavior where x0=0 and hx=0.2 gives centers like 5.3.
    i = pos[1]
    j = pos[2]
    k = pos[3]

    cx = x0g[1] + (Float64(i) - 0.5) * hxg[1]
    cy = x0g[2] + (Float64(j) - 0.5) * hxg[2]
    cθ = x0g[3] + (Float64(k) - 0.5) * hxg[3]

    return (cx, cy, cθ)
end

"""Fast (cached) 2D center of an abstract state."""
function center2d(s::Int)
    if haskey(_CENTER2D_CACHE, s)
        return _CENTER2D_CACHE[s]
    end

    cx, cy, _ = _center3d_from_grid(abstract_system, s)
    xy = (Float64(cx), Float64(cy))
    _CENTER2D_CACHE[s] = xy
    return xy
end


if !isdefined(Main, :_CENTER2D_CACHE)
    const _CENTER2D_CACHE = Dict{Int, Tuple{Float64, Float64}}()
end

