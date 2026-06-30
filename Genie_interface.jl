

############################################################
# Genie web interface for the interactive symbolic planner
#
# Usage from the Julia REPL:
#   include("Interactive_GPT.jl")
#   include("Genie_interface.jl")
#   start_genie_interface()
#
# Then open:
#   http://127.0.0.1:8000
############################################################

using Genie
using Genie.Router
using Genie.Renderer.Json
using Genie.Requests
using JSON3

const GENIE_DEFAULT_HOST = "127.0.0.1"
const GENIE_DEFAULT_PORT = 8000
const GENIE_SUBMIT_LOCK = ReentrantLock()
const GENIE_PENDING_BUILDS = IdDict{Any,Any}()
const GENIE_PENDING_NL = IdDict{Any,String}()
using Genie.Renderer.Html
"""Return a named object from Main if it exists, otherwise `nothing`."""
function _main_get(name::Symbol)
    return isdefined(Main, name) ? getfield(Main, name) : nothing
end


"""Try to read a field from an object, returning `nothing` if unavailable."""
function _maybe_getfield(obj, field::Symbol)
    obj === nothing && return nothing
    return hasfield(typeof(obj), field) ? getfield(obj, field) : nothing
end

"""Convert a workspace rectangle in [0,10]x[0,10] to CSS absolute-position percentages.
The workspace origin is bottom-left, while CSS uses top-left, so the y-axis is inverted.
"""
function _css_rect_from_workspace(xmin, xmax, ymin, ymax)
    left = 10.0 * Float64(xmin)
    top = 100.0 - 10.0 * Float64(ymax)
    width = 10.0 * (Float64(xmax) - Float64(xmin))
    height = 10.0 * (Float64(ymax) - Float64(ymin))
    return "left: $(round(left; digits=2))%; top: $(round(top; digits=2))%; width: $(round(width; digits=2))%; height: $(round(height; digits=2))%;"
end

"""Convert a workspace point in [0,10]x[0,10] to CSS absolute-position percentages."""

function _css_point_from_workspace(x, y)
    left = 10.0 * Float64(x)
    top = 100.0 - 10.0 * Float64(y)
    return "left: $(round(left; digits=2))%; top: $(round(top; digits=2))%;"
end

"""Return the quotient-state id associated with a workspace point."""
function _quotient_state_from_xy(x, y)
    x = Float64(x)
    y = Float64(y)

    if 6.0 <= x <= 9.0 && 7.5 <= y <= 9.5
        if 4.0 <= x <= 7.5 && 7.0 <= y <= 8.0
            return "q3"   # {blue,brown}
        end
        return "q2"       # {blue}
    elseif 4.0 <= x <= 7.5 && 7.0 <= y <= 8.0
        return "q4"       # {brown}
    elseif 7.0 <= x <= 8.5 && 3.5 <= y <= 6.5
        return "q5"       # {purple}
    elseif 5.0 <= x <= 7.5 && 1.0 <= y <= 2.5
        return "q6"       # {green}
    elseif 0.5 <= x <= 2.5 && 1.0 <= y <= 3.5
        return "q7"       # {yellow}
    end

    if isdefined(Main, :x1_lb) && isdefined(Main, :x1_ub) && isdefined(Main, :x2_lb) && isdefined(Main, :x2_ub)
        x1_lb = getfield(Main, :x1_lb)
        x1_ub = getfield(Main, :x1_ub)
        x2_lb = getfield(Main, :x2_lb)
        x2_ub = getfield(Main, :x2_ub)
        for i in eachindex(x1_lb)
            if Float64(x1_lb[i]) <= x <= Float64(x1_ub[i]) && Float64(x2_lb[i]) <= y <= Float64(x2_ub[i])
                return "q1"   # {obstacle}
            end
        end
    end

    return "q0"           # free unlabeled space
end

"""Generate the workspace map HTML using the same region/obstacle geometry as the planner."""
function genie_workspace_map_html()
    region_html = String[]

    # AP regions from the current benchmark environment.
    push!(region_html, "<div class=\"region r-blue\" style=\"$(_css_rect_from_workspace(6.0, 9.0, 7.5, 9.5))\"></div>")
    push!(region_html, "<div class=\"region r-green\" style=\"$(_css_rect_from_workspace(5.0, 7.5, 1.0, 2.5))\"></div>")
    push!(region_html, "<div class=\"region r-purple\" style=\"$(_css_rect_from_workspace(7.0, 8.5, 3.5, 6.5))\"></div>")
    push!(region_html, "<div class=\"region r-brown\" style=\"$(_css_rect_from_workspace(4.0, 7.5, 7.0, 8.0))\"></div>")
    push!(region_html, "<div class=\"region r-yellow\" style=\"$(_css_rect_from_workspace(0.5, 2.5, 1.0, 3.5))\"></div>")

    # Obstacles from the planner globals when available.
    if isdefined(Main, :x1_lb) && isdefined(Main, :x1_ub) && isdefined(Main, :x2_lb) && isdefined(Main, :x2_ub)
        x1_lb = getfield(Main, :x1_lb)
        x1_ub = getfield(Main, :x1_ub)
        x2_lb = getfield(Main, :x2_lb)
        x2_ub = getfield(Main, :x2_ub)
        for i in eachindex(x1_lb)
            push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(x1_lb[i], x1_ub[i], x2_lb[i], x2_ub[i]))\"></div>")
        end
    else
        # Fallback obstacle geometry matching the current workspace draft.
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(0.5, 2.5, 6.5, 9.0))\"></div>")
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(0.5, 2.5, 5.0, 5.5))\"></div>")
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(3.5, 4.0, 3.5, 6.0))\"></div>")
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(5.0, 5.5, 3.5, 6.0))\"></div>")
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(3.5, 5.5, 3.5, 4.0))\"></div>")
        push!(region_html, "<div class=\"region obs\" style=\"$(_css_rect_from_workspace(3.5, 4.0, 1.0, 3.0))\"></div>")
    end

    vehicle_x, vehicle_y = 4.8, 5.0
    runtime = _main_get(:RUNTIME)
    if runtime !== nothing && hasfield(typeof(runtime), :current_state)
        vehicle_x = runtime.current_state[1]
        vehicle_y = runtime.current_state[2]
    end

    push!(region_html, "<div class=\"trajectory-layer\"><svg viewBox=\"0 0 100 100\" preserveAspectRatio=\"none\"><polyline id=\"trajectory-polyline\" points=\"\"/></svg></div>")
    push!(region_html, "<div id=\"vehicle\" class=\"vehicle\" style=\"$(_css_point_from_workspace(vehicle_x, vehicle_y))\"></div>")

    return join(region_html, "\n        ")
end

"""Static quotient-transition-system schematic for the current workspace partition.
This groups abstract states by their visible atomic-proposition label sets.
"""
function genie_quotient_svg_html()
    return """
<svg class=\"quotient-svg\" viewBox=\"0 0 760 430\" xmlns=\"http://www.w3.org/2000/svg\">
  <defs>
    <marker id=\"q-arrow\" markerWidth=\"14\" markerHeight=\"14\" refX=\"11\" refY=\"7\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">
      <path d=\"M0,0 L0,14 L12,7 z\" fill=\"#0f172a\" />
    </marker>
    <marker id=\"q-arrow-start\" markerWidth=\"14\" markerHeight=\"14\" refX=\"1\" refY=\"7\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">
      <path d=\"M12,0 L12,14 L0,7 z\" fill=\"#0f172a\" />
    </marker>
  </defs>

  <style>
    .qedge {
      fill: none;
      stroke: #0f172a;
      stroke-width: 2.8;
      stroke-linecap: round;
      stroke-linejoin: round;
      marker-end: url(#q-arrow);
    }
    .qedge-bidir {
      fill: none;
      stroke: #0f172a;
      stroke-width: 2.8;
      stroke-linecap: round;
      stroke-linejoin: round;
      marker-start: url(#q-arrow-start);
      marker-end: url(#q-arrow);
    }
    .qloop {
      fill: none;
      stroke: #0f172a;
      stroke-width: 2.8;
      stroke-linecap: round;
      stroke-linejoin: round;
      marker-end: url(#q-arrow);
    }
    .qnode {
      fill: #ffffff;
      stroke: #0f172a;
      stroke-width: 2.8;
    }
    .qactive {
      fill: #dbeafe;
      stroke: #2563eb;
      stroke-width: 4.0;
    }
    .qlabel {
      font: 700 20px Inter, sans-serif;
      fill: #111827;
      text-anchor: middle;
      dominant-baseline: middle;
    }
    .aplabel {
      font: 16px Inter, sans-serif;
      fill: #111827;
      text-anchor: middle;
    }
  </style>

  <!-- Main outgoing and return edges as bidirectional edges. -->
  <path class=\"qedge-bidir\" d=\"M350,192 C350,135 360,105 370,105\" />
  <path class=\"qedge-bidir\" d=\"M412,78 C465,62 515,78 542,105\" />
  <path class=\"qedge-bidir\" d=\"M370,205 C440,145 500,120 535,120\" />
  <path class=\"qedge-bidir\" d=\"M375,225 C470,220 560,220 620,220\" />
  <path class=\"qedge-bidir\" d=\"M365,245 C455,295 535,325 585,325\" />
  <path class=\"qedge-bidir\" d=\"M342,260 C340,305 350,330 360,333\" />
  <path class=\"qedge-bidir\" d=\"M122,220 C190,220 245,220 305,220\" />
  <path class=\"qedge-bidir\" d=\"M315,252 C260,295 220,330 180,360\" />

  <!-- Brown/overlap local connections as bidirectional edge. -->
  <path class=\"qedge-bidir\" d=\"M595,135 C625,155 640,178 645,190\" />

  <!-- Self loops, drawn outside the nodes so they are clearly visible. -->
  <path class=\"qloop\" d=\"M70,185 C25,155 25,285 70,255\" />
  <path class=\"qloop\" d=\"M350,48 C315,0 445,0 410,48\" />
  <path class=\"qloop\" d=\"M585,92 C620,58 680,105 622,138\" />
  <path class=\"qloop\" d=\"M682,195 C740,165 740,275 682,245\" />
  <path class=\"qloop\" d=\"M642,310 C705,300 705,390 642,360\" />
  <path class=\"qloop\" d=\"M350,392 C315,445 430,445 395,392\" />
  <path class=\"qloop\" d=\"M126,397 C70,430 100,500 168,405\" />
  <path class=\"qloop\" d=\"M320,190 C260,145 260,305 320,260\" />

  <!-- Nodes. -->
  <circle id="quotient-q0" class="qnode qactive" cx="340" cy="225" r="34" />
  <circle id="quotient-q1" class="qnode" cx="90" cy="220" r="32" />
  <circle id="quotient-q2" class="qnode" cx="380" cy="75" r="32" />
  <circle id="quotient-q3" class="qnode" cx="570" cy="125" r="32" />
  <circle id="quotient-q4" class="qnode" cx="655" cy="220" r="32" />
  <circle id="quotient-q5" class="qnode" cx="615" cy="335" r="32" />
  <circle id="quotient-q6" class="qnode" cx="365" cy="365" r="32" />
  <circle id="quotient-q7" class="qnode" cx="150" cy="380" r="32" />

  <!-- State labels. -->
  <text class=\"qlabel\" x=\"340\" y=\"225\">q₀</text>
  <text class=\"qlabel\" x=\"90\" y=\"220\">q₁</text>
  <text class=\"qlabel\" x=\"380\" y=\"75\">q₂</text>
  <text class=\"qlabel\" x=\"570\" y=\"125\">q₃</text>
  <text class=\"qlabel\" x=\"655\" y=\"220\">q₄</text>
  <text class=\"qlabel\" x=\"615\" y=\"335\">q₅</text>
  <text class=\"qlabel\" x=\"365\" y=\"365\">q₆</text>
  <text class=\"qlabel\" x=\"150\" y=\"380\">q₇</text>

  <!-- AP labels. -->
  <text class="aplabel" x="88" y="156">{obstacle}</text>
  <text class="aplabel" x="382" y="0">{blue}</text>
  <text class="aplabel" x="570" y="40">{blue,brown}</text>
  <text class="aplabel" x="708" y="172">{brown}</text>
  <text class="aplabel" x="690" y="384">{purple}</text>
  <text class="aplabel" x="375" y="455">{green}</text>
  <text class="aplabel" x="150" y="340">{yellow}</text>
</svg>
"""
end



"""Convert different state representations to a browser plottable [x1,x2] point."""
function _xy_point_from_state_like(z)
    z === nothing && return nothing

    # Continuous state, e.g. SVector(x1,x2,theta), Vector, Tuple.
    try
        if length(z) >= 2
            return [Float64(z[1]), Float64(z[2])]
        end
    catch
    end

    # Abstract state index: use existing center2d if available.
    center2d = _main_get(:center2d)
    if center2d !== nothing
        try
            x, y = center2d(z)
            return [Float64(x), Float64(y)]
        catch
        end
    end

    return nothing
end


"""Return the planner integration step in seconds for browser animation."""
function _planner_dt_seconds()
    return isdefined(Main, :DT) ? Float64(getfield(Main, :DT)) : 0.30
end

"""Advance the persistent runtime to the end of the latest solved trajectory.
This makes the next submitted requirement start from the vehicle's current location,
not from the original initial condition.
"""
function _advance_runtime_to_latest_plan_end!(runtime)
    runtime === nothing && return nothing
    isdefined(Main, :_RUNTIME_PLAN_XS) || return nothing

    xs_dict = getfield(Main, :_RUNTIME_PLAN_XS)
    haskey(xs_dict, runtime) || return nothing

    xs = xs_dict[runtime]
    isempty(xs) && return nothing

    final_abs = xs[end]

    try
        runtime.current_abs = final_abs
    catch
    end

    center2d = _main_get(:center2d)
    if center2d !== nothing && hasfield(typeof(runtime), :current_state)
        try
            x, y = center2d(final_abs)
            old_state = runtime.current_state
            θ = length(old_state) >= 3 ? Float64(old_state[3]) : 0.0
            runtime.current_state = typeof(old_state)(Float64(x), Float64(y), θ)
        catch
            try
                x, y = center2d(final_abs)
                runtime.current_state = [Float64(x), Float64(y), 0.0]
            catch
            end
        end
    end

    return final_abs
end

"""Reset the runtime start state to the vehicle position currently displayed in the browser.
This is used when the user interrupts an animation with a new requirement: the new plan
should start from the displayed stop position, not from the previous plan destination.
"""
function _set_runtime_current_state_from_xy!(runtime, xy)
    runtime === nothing && return nothing
    xy === nothing && return nothing

    x = nothing
    y = nothing
    try
        if length(xy) >= 2
            x = Float64(xy[1])
            y = Float64(xy[2])
        end
    catch
        return nothing
    end
    x === nothing && return nothing
    y === nothing && return nothing

    if hasfield(typeof(runtime), :current_state)
        try
            old_state = runtime.current_state
            θ = length(old_state) >= 3 ? Float64(old_state[3]) : 0.0
            try
                runtime.current_state = typeof(old_state)(x, y, θ)
            catch
                runtime.current_state = [x, y, θ]
            end
        catch
        end
    end

    # Best effort: update the abstract current state as the nearest state center.
    # This keeps symbolic planning aligned with the browser stop position when
    # center2d and the transition-system state count are available.
    center2d = _main_get(:center2d)
    backend = _backend_from_runtime(runtime)
    if center2d !== nothing && backend !== nothing && hasfield(typeof(backend), :ts)
        try
            ts = backend.ts
            nstates = 0

            L2A = _main_get(:L2A)
            if L2A !== nothing
                try
                    nstates = Int(L2A.nstates(ts))
                catch
                end
            end

            if nstates <= 0 && hasfield(typeof(ts), :n)
                nstates = Int(getfield(ts, :n))
            elseif nstates <= 0 && hasfield(typeof(ts), :nx)
                nstates = Int(getfield(ts, :nx))
            elseif nstates <= 0 && hasfield(typeof(ts), :X)
                nstates = length(getfield(ts, :X))
            elseif nstates <= 0 && hasfield(typeof(ts), :Adj)
                nstates = length(getfield(ts, :Adj))
            elseif nstates <= 0 && hasfield(typeof(ts), :adj)
                nstates = length(getfield(ts, :adj))
            end

            nstates <= 0 && return nothing

            best_abs = nothing
            best_dist = Inf
            for abs_state in 1:nstates
                try
                    cx, cy = center2d(abs_state)
                    d = (Float64(cx) - x)^2 + (Float64(cy) - y)^2
                    if d < best_dist
                        best_dist = d
                        best_abs = abs_state
                    end
                catch
                end
            end

            if best_abs !== nothing
                if hasfield(typeof(runtime), :current_abs)
                    try
                        runtime.current_abs = best_abs
                    catch
                    end
                end
                if hasfield(typeof(runtime), :trajectory_abs)
                    try
                        runtime.trajectory_abs = [best_abs]
                    catch
                    end
                end
                if hasfield(typeof(runtime), :trajectory_xy)
                    try
                        runtime.trajectory_xy = [(x, y)]
                    catch
                    end
                end
                @info "Runtime current state reset from browser position" x=x y=y current_abs=best_abs nearest_center_distance=sqrt(best_dist)
                return best_abs
            end
        catch
        end
    end

    return nothing
end

"""Read the actual trajectory produced by the symbolic controller/runtime.
This intentionally does not fabricate a path from the formula. It only returns data already computed by the planner.
"""
function _actual_runtime_trajectory(runtime)
    runtime === nothing && return Vector{Vector{Float64}}()

    # 1) Preferred: dictionaries used by the GLMakie GUI runtime.
    for dict_name in (:_RUNTIME_PLAN_XS, :_RUNTIME_PLAN_STATES, :_RUNTIME_PLAN_ABS, :_RUNTIME_TRAJECTORY)
        if isdefined(Main, dict_name)
            d = getfield(Main, dict_name)
            try
                if haskey(d, runtime)
                    raw = d[runtime]
                    pts = Vector{Vector{Float64}}()
                    for z in raw
                        p = _xy_point_from_state_like(z)
                        p === nothing || push!(pts, p)
                    end
                    length(pts) >= 1 && return pts
                end
            catch
            end
        end
    end

    # 2) Common runtime fields, depending on local script version.
    for field in (:xs, :plan_xs, :current_xs, :trajectory, :current_trajectory, :executed_states, :path, :plan)
        if hasfield(typeof(runtime), field)
            try
                raw = getfield(runtime, field)
                pts = Vector{Vector{Float64}}()
                for z in raw
                    p = _xy_point_from_state_like(z)
                    p === nothing || push!(pts, p)
                end
                length(pts) >= 1 && return pts
            catch
            end
        end
    end

    # 3) Fallback to current state only; no fake future path.
    if hasfield(typeof(runtime), :current_state)
        p = _xy_point_from_state_like(runtime.current_state)
        p === nothing || return [p]
    end

    return Vector{Vector{Float64}}()
end

"""Try to run the same non-GUI build/planning pipeline used by the GLMakie frontend."""
function _run_existing_symbolic_planner!(runtime, nl::AbstractString)
    runtime === nothing && error("RUNTIME is not available. Include Interface_GENIE.jl or Interactive_GPT.jl first.")

    # If a project-specific non-GUI entry point exists, use it first.
    for fname in (:run_buchi_from_user_nl!, :run_buchi_from_nl!, :run_from_user_nl!)
        f = _main_get(fname)
        f === nothing && continue
        try
            return f(String(nl))
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            @warn "Planner entry point failed; trying another signature" function_name=fname error=msg
            try
                return f(runtime, String(nl))
            catch err2
                msg2 = sprint(showerror, err2)
                if occursin("quotient", lowercase(msg2)) ||
                   occursin("not compatible", lowercase(msg2)) ||
                   occursin("infeasible", lowercase(msg2)) ||
                   occursin("no_admissible_control", lowercase(msg2))
                    rethrow(err2)
                end
            end
        end
    end

    # Otherwise reproduce the GLMakie submit path using the existing low-level helpers.
    build_fns = (:build_buchi_from_user_nl, :build_buchi_from_user_nl!, :build_buchi_from_nl, :build_buchi_from_nl!)
    build_result = nothing

    backend = _maybe_getfield(runtime, :backend)
    backend === nothing && (backend = _main_get(:BACKEND))
    backend === nothing && (backend = _main_get(:backend))

    for fname in build_fns
        f = _main_get(fname)
        f === nothing && continue

        # Preferred current signature: build_buchi_from_nl(backend, nl)
        if backend !== nothing
            try
                build_result = f(backend, String(nl))
                break
            catch err
                msg = sprint(showerror, err)
                if occursin("quotient", lowercase(msg)) ||
                   occursin("not compatible", lowercase(msg)) ||
                   occursin("infeasible", lowercase(msg)) ||
                   occursin("no_admissible_control", lowercase(msg))
                    rethrow(err)
                end
                @warn "Büchi build entry point failed; trying another signature" function_name=fname error=msg
            end
        end

        # Fallback signatures from older versions.
        try
            build_result = f(runtime, String(nl))
            break
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            @warn "Büchi build entry point failed; trying another signature" function_name=fname error=msg
        end

        try
            build_result = f(String(nl))
            break
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            @warn "Büchi build entry point failed; trying next build function" function_name=fname error=msg
        end
    end

    if build_result !== nothing
        reset_runtime_plan! = _main_get(:reset_runtime_plan!)
        install_buchi_build! = _main_get(:install_buchi_build!)
        plan_from_runtime! = _main_get(:plan_from_runtime!)

        reset_runtime_plan! === nothing || reset_runtime_plan!(runtime)
        install_buchi_build! === nothing || install_buchi_build!(runtime, build_result)

        if plan_from_runtime! !== nothing
            return plan_from_runtime!(runtime)
        end
        return build_result
    end

    # Last resort: translate only. The returned trajectory will still be current-state only.
    return nothing
end

function _backend_from_runtime(runtime)
    backend = _maybe_getfield(runtime, :backend)
    backend === nothing && (backend = _main_get(:BACKEND))
    backend === nothing && (backend = _main_get(:backend))
    return backend
end

function _build_buchi_for_validation(runtime, nl::AbstractString)
    backend = _backend_from_runtime(runtime)
    backend === nothing && error("Planner backend is not available. Build the symbolic backend before starting Genie.")

    build_fns = (:build_buchi_from_user_nl, :build_buchi_from_user_nl!, :build_buchi_from_nl, :build_buchi_from_nl!)
    last_error = nothing

    for fname in build_fns
        f = _main_get(fname)
        f === nothing && continue

        try
            return f(backend, String(nl))
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            last_error = err
            @warn "Validation build entry point failed; trying another signature" function_name=fname error=msg
        end

        try
            return f(runtime, String(nl))
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            last_error = err
            @warn "Validation build entry point failed; trying another signature" function_name=fname error=msg
        end

        try
            return f(String(nl))
        catch err
            msg = sprint(showerror, err)
            if occursin("quotient", lowercase(msg)) ||
               occursin("not compatible", lowercase(msg)) ||
               occursin("infeasible", lowercase(msg)) ||
               occursin("no_admissible_control", lowercase(msg))
                rethrow(err)
            end
            last_error = err
            @warn "Validation build entry point failed; trying next build function" function_name=fname error=msg
        end
    end

    last_error === nothing || throw(last_error)
    error("No Büchi build function is available for validation.")
end

function _quotient_check_payload()
    quotient_message = ""
    quotient_product_states = nothing
    quotient_states_count = nothing
    if isdefined(Main, :LAST_QUOTIENT_CHECK)
        qc = Main.LAST_QUOTIENT_CHECK
        try
            quotient_message = String(getfield(qc, :message))
            quotient_product_states = getfield(qc, :nproduct_states_reached)
            quotient_states_count = getfield(qc, :nquotient_states)
        catch
            quotient_message = string(qc)
        end
    end
    return quotient_message, quotient_product_states, quotient_states_count
end

function genie_prepare_validation(nl_text::AbstractString, current_state_xy = nothing)
    nl = strip(String(nl_text))
    isempty(nl) && error("Empty natural-language requirement.")

    runtime = _main_get(:RUNTIME)
    runtime === nothing && error("RUNTIME is not available. Include the interactive planner file before starting Genie.")

    lock(GENIE_SUBMIT_LOCK)
    try
        begin_replan = _main_get(:begin_replan!)
        if begin_replan !== nothing
            begin_replan(runtime, nl)
        else
            token = _maybe_getfield(runtime, :plan_token)
            if token !== nothing
                try
                    setfield!(runtime, :plan_token, token + 1)
                catch
                end
            end
        end
        delete!(GENIE_PENDING_BUILDS, runtime)
        delete!(GENIE_PENDING_NL, runtime)

        reset_runtime_plan = _main_get(:reset_runtime_plan!)
        reset_runtime_plan === nothing || reset_runtime_plan(runtime)

        stopped_abs = _set_runtime_current_state_from_xy!(runtime, current_state_xy)
        if stopped_abs !== nothing
            @info "New requirement starts from interrupted browser position" current_state_xy=current_state_xy current_abs=stopped_abs
        elseif current_state_xy !== nothing
            @info "New requirement received browser stop position; updated continuous runtime state only" current_state_xy=current_state_xy
        end

        build_result = _build_buchi_for_validation(runtime, nl)
        GENIE_PENDING_BUILDS[runtime] = build_result
        GENIE_PENDING_NL[runtime] = nl

        formula = hasfield(typeof(build_result), :formula) ? String(getfield(build_result, :formula)) : ""
        backtranslation = hasfield(typeof(build_result), :backtranslation) ? String(getfield(build_result, :backtranslation)) : ""
        validation_message = hasfield(typeof(build_result), :validation_message) ? String(getfield(build_result, :validation_message)) : ""
        validation_message = replace(validation_message, r"^\s*\*\*?Validation Message:?\*\*?\s*" => "")
        validation_message = replace(validation_message, r"^\s*Validation Message:?\s*" => "")
        isempty(validation_message) && (validation_message = "I understood your request as: $(isempty(backtranslation) ? formula : backtranslation). Is this correct?")

        quotient_message, quotient_product_states, quotient_states_count = _quotient_check_payload()
        status_text = validation_message
        if !isempty(quotient_message) && strip(quotient_message) != "Requirement passed the coarse AP-label quotient compatibility check."
            status_text *= "\n\n" * quotient_message
            if quotient_states_count !== nothing && quotient_product_states !== nothing
                status_text *= "\nQuotient states: " * string(quotient_states_count) *
                               ", reachable product states: " * string(quotient_product_states)
            end
        end

        return Dict(
            "ok" => true,
            "needs_validation" => true,
            "plan_token" => _maybe_getfield(runtime, :plan_token),
            "nl" => nl,
            "formula" => formula,
            "backtranslation" => backtranslation,
            "validation_message" => validation_message,
            "quotient_message" => quotient_message,
            "quotient_state_count" => quotient_states_count,
            "quotient_product_state_count" => quotient_product_states,
            "status" => status_text,
        )
    finally
        unlock(GENIE_SUBMIT_LOCK)
    end
end

function genie_confirm_validation_and_plan()
    runtime = _main_get(:RUNTIME)
    runtime === nothing && error("RUNTIME is not available. Include the interactive planner file before starting Genie.")

    lock(GENIE_SUBMIT_LOCK)
    try
        haskey(GENIE_PENDING_BUILDS, runtime) || error("There is no pending validated requirement. Submit a requirement first.")
        build_result = GENIE_PENDING_BUILDS[runtime]
        nl = get(GENIE_PENDING_NL, runtime, "")
        token_before_plan = _maybe_getfield(runtime, :plan_token)

        reset_runtime_plan = _main_get(:reset_runtime_plan!)
        install_buchi_build = _main_get(:install_buchi_build!)
        plan_from_runtime = _main_get(:plan_from_runtime!)

        reset_runtime_plan === nothing || reset_runtime_plan(runtime)
        install_buchi_build === nothing && error("install_buchi_build! is not available.")
        plan_from_runtime === nothing && error("plan_from_runtime! is not available.")

        install_buchi_build(runtime, build_result)
        planning_status = plan_from_runtime(runtime)
        token_after_plan = _maybe_getfield(runtime, :plan_token)
        if token_before_plan !== nothing && token_after_plan !== nothing && token_before_plan != token_after_plan
            error("A newer requirement was submitted while this plan was being computed. Ignoring the stale plan.")
        end

        delete!(GENIE_PENDING_BUILDS, runtime)
        delete!(GENIE_PENDING_NL, runtime)

        formula = hasfield(typeof(build_result), :formula) ? String(getfield(build_result, :formula)) : ""
        backtranslation = hasfield(typeof(build_result), :backtranslation) ? String(getfield(build_result, :backtranslation)) : ""
        validation_message = hasfield(typeof(build_result), :validation_message) ? String(getfield(build_result, :validation_message)) : ""
        validation_message = replace(validation_message, r"^\s*\*\*?Validation Message:?\*\*?\s*" => "")
        validation_message = replace(validation_message, r"^\s*Validation Message:?\s*" => "")

        quotient_message, quotient_product_states, quotient_states_count = _quotient_check_payload()

        trajectory = _actual_runtime_trajectory(runtime)
        # Do not advance the persistent runtime to the planned destination here.
        # The browser animation is the execution clock. If the user interrupts,
        # /api/submit sends the displayed vehicle position back and the runtime is
        # reset from that actual stop position.
        quotient_states = [_quotient_state_from_xy(p[1], p[2]) for p in trajectory]

        status_text = "Validation accepted. Solved with the existing abstraction-based symbolic controller. Status: " * string(planning_status)
        if !isempty(quotient_message) && strip(quotient_message) != "Requirement passed the coarse AP-label quotient compatibility check."
            status_text *= "\n" * quotient_message
            if quotient_states_count !== nothing && quotient_product_states !== nothing
                status_text *= "\nQuotient states: " * string(quotient_states_count) *
                               ", reachable product states: " * string(quotient_product_states)
            end
        end

        return Dict(
            "ok" => true,
            "needs_validation" => false,
            "plan_token" => _maybe_getfield(runtime, :plan_token),
            "nl" => nl,
            "formula" => formula,
            "backtranslation" => backtranslation,
            "validation_message" => validation_message,
            "trajectory" => trajectory,
            "dt_seconds" => _planner_dt_seconds(),
            "quotient_states" => quotient_states,
            "quotient_message" => quotient_message,
            "quotient_state_count" => quotient_states_count,
            "quotient_product_state_count" => quotient_product_states,
            "status" => status_text,
        )
    finally
        unlock(GENIE_SUBMIT_LOCK)
    end
end

function genie_reject_validation()
    runtime = _main_get(:RUNTIME)
    runtime === nothing && error("RUNTIME is not available. Include the interactive planner file before starting Genie.")

    lock(GENIE_SUBMIT_LOCK)
    try
        begin_replan = _main_get(:begin_replan!)
        if begin_replan !== nothing
            begin_replan(runtime, "")
        end
        delete!(GENIE_PENDING_BUILDS, runtime)
        delete!(GENIE_PENDING_NL, runtime)
        return Dict(
            "ok" => true,
            "needs_validation" => false,
            "status" => "Validation rejected. Please rephrase the requirement and submit again.",
        )
    finally
        unlock(GENIE_SUBMIT_LOCK)
    end
end

function genie_translate_and_plan(nl_text::AbstractString)
    nl = strip(String(nl_text))
    isempty(nl) && error("Empty natural-language requirement.")

    runtime = _main_get(:RUNTIME)
    runtime === nothing && error("RUNTIME is not available. Include the interactive planner file before starting Genie.")

    lock(GENIE_SUBMIT_LOCK)
    try
        # Run the real symbolic-control pipeline. This should update the runtime with
        # the Spot automaton, ABCD/symbolic plan, quotient check, and controller trajectory.
        planning_status = _run_existing_symbolic_planner!(runtime, nl)

        # Some planner entry points update RUNTIME in-place but return `nothing`.
        # Treat this as valid if the runtime/global planner dictionaries now contain
        # a Büchi build and a synthesized/simulated plan.
        has_buchi_build = _maybe_getfield(runtime, :current_aut) !== nothing &&
                          _maybe_getfield(runtime, :current_labels_ts) !== nothing
        has_plan_status = isdefined(Main, :_RUNTIME_PLAN_STATUS) && haskey(Main._RUNTIME_PLAN_STATUS, runtime)
        has_plan_path = isdefined(Main, :_RUNTIME_PLAN_XS) && haskey(Main._RUNTIME_PLAN_XS, runtime) &&
                        !isempty(get(Main._RUNTIME_PLAN_XS, runtime, Int[]))

        if planning_status === nothing && !has_plan_status && !has_plan_path && has_buchi_build
            plan_from_runtime = _main_get(:plan_from_runtime!)
            if plan_from_runtime !== nothing
                planning_status = plan_from_runtime(runtime)
                has_plan_status = isdefined(Main, :_RUNTIME_PLAN_STATUS) && haskey(Main._RUNTIME_PLAN_STATUS, runtime)
                has_plan_path = isdefined(Main, :_RUNTIME_PLAN_XS) && haskey(Main._RUNTIME_PLAN_XS, runtime) &&
                                !isempty(get(Main._RUNTIME_PLAN_XS, runtime, Int[]))
            end
        end

        if planning_status === nothing && !(has_buchi_build && (has_plan_status || has_plan_path))
            # Some entry points build LAST_BUCHI_BUILD globally but do not install
            # it into RUNTIME. Install it here and synthesize once.
            if isdefined(Main, :LAST_BUCHI_BUILD)
                install_buchi_build = _main_get(:install_buchi_build!)
                reset_runtime_plan = _main_get(:reset_runtime_plan!)
                plan_from_runtime = _main_get(:plan_from_runtime!)

                if reset_runtime_plan !== nothing
                    reset_runtime_plan(runtime)
                end

                if install_buchi_build !== nothing
                    install_buchi_build(runtime, Main.LAST_BUCHI_BUILD)
                    has_buchi_build = _maybe_getfield(runtime, :current_aut) !== nothing &&
                                      _maybe_getfield(runtime, :current_labels_ts) !== nothing
                end

                if plan_from_runtime !== nothing && has_buchi_build
                    planning_status = plan_from_runtime(runtime)
                    has_plan_status = isdefined(Main, :_RUNTIME_PLAN_STATUS) && haskey(Main._RUNTIME_PLAN_STATUS, runtime)
                    has_plan_path = isdefined(Main, :_RUNTIME_PLAN_XS) && haskey(Main._RUNTIME_PLAN_XS, runtime) &&
                                    !isempty(get(Main._RUNTIME_PLAN_XS, runtime, Int[]))
                end
            end
        end

        if planning_status === nothing && !(has_buchi_build && (has_plan_status || has_plan_path))
            error("No symbolic plan was produced. The request may be incompatible with the system capabilities or the planner did not return a result.")
        end

        if planning_status === nothing && has_plan_status
            planning_status = get(Main._RUNTIME_PLAN_STATUS, runtime, nothing)
        end

        quotient_message = ""
        quotient_product_states = nothing
        quotient_states_count = nothing
        if isdefined(Main, :LAST_QUOTIENT_CHECK)
            qc = Main.LAST_QUOTIENT_CHECK
            try
                quotient_message = String(getfield(qc, :message))
                quotient_product_states = getfield(qc, :nproduct_states_reached)
                quotient_states_count = getfield(qc, :nquotient_states)
            catch
                quotient_message = string(qc)
            end
        end

        formula = ""
        if hasfield(typeof(runtime), :current_formula)
            formula = String(runtime.current_formula)
        end

        # If the runtime did not store the formula, compute it only for display.
        if isempty(formula)
            nl_to_spot_formula = _main_get(:nl_to_spot_formula)
            backend = _maybe_getfield(runtime, :backend)
            backend === nothing && (backend = _main_get(:BACKEND))
            backend === nothing && (backend = _main_get(:backend))
            if nl_to_spot_formula !== nothing && backend !== nothing
                formula = String(nl_to_spot_formula(backend, nl))
            end
        end
        if !has_plan_path && isdefined(Main, :_RUNTIME_PLAN_XS) && !haskey(Main._RUNTIME_PLAN_XS, runtime)
            error("Planner completed but no runtime trajectory was stored for the interface.")
        end
        trajectory = _actual_runtime_trajectory(runtime)
        # Legacy direct planning path: do not force the runtime to the destination.
        # The validation interface updates the runtime from the displayed browser position.
        quotient_states = [_quotient_state_from_xy(p[1], p[2]) for p in trajectory]

        status_text = "Solved with the existing abstraction-based symbolic controller. Status: " *
                      string(planning_status === nothing ? "updated runtime in-place" : planning_status)
        if !isempty(quotient_message)
            status_text *= "\n" * quotient_message
            if quotient_states_count !== nothing && quotient_product_states !== nothing
                status_text *= "\nQuotient states: " * string(quotient_states_count) *
                               ", reachable product states: " * string(quotient_product_states)
            end
        end

        return Dict(
            "ok" => true,
            "nl" => nl,
            "formula" => formula,
            "trajectory" => trajectory,
            "dt_seconds" => _planner_dt_seconds(),
            "quotient_states" => quotient_states,
            "quotient_message" => quotient_message,
            "quotient_state_count" => quotient_states_count,
            "quotient_product_state_count" => quotient_product_states,
            "status" => status_text,
        )
    finally
        unlock(GENIE_SUBMIT_LOCK)
    end
end

"""HTML dashboard. Kept in Julia for easy single-file experimentation."""
function genie_index_html()
    return """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Interactive Symbolic Planner</title>
  <style>
    :root {
      --bg: #64727d;
      --card: #ffffff;
      --ink: #111827;
      --muted: #6b7280;
      --line: #111827;
      --blue: #2563eb;
      --green: #16a34a;
      --purple: #7c3aed;
      --red: #dc2626;
      --yellow: #facc15;
      --brown: #7c2d12;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--ink);
    }
    .app {
      width: min(1440px, 96vw);
      height: min(900px, 94vh);
      margin: 3vh auto;
      display: grid;
      grid-template-columns: 1.45fr 0.95fr 0.52fr;
      grid-template-rows: 0.9fr 1.78fr 0.82fr 0.65fr;
      gap: 14px;
    }
    .card {
      background: var(--card);
      border: 2px solid var(--line);
      border-radius: 14px;
      box-shadow: 0 14px 30px rgba(0, 0, 0, 0.14);
      padding: 16px;
      overflow: hidden;
    }
    .title {
      font-size: 20px;
      font-weight: 800;
      text-align: center;
      margin-bottom: 12px;
    }
    .workspace { grid-row: 1 / span 3; grid-column: 1; padding: 14px; }
    .formula { grid-row: 1; grid-column: 2; }
    .automaton { grid-row: 2; grid-column: 2; }
    .ap { grid-row: 3; grid-column: 2; margin-top: 4px; }
    .ap .title { margin-bottom: 8px; }
    .info { grid-row: 2; grid-column: 2; }
    .input { grid-row: 4; grid-column: 1; display: flex; align-items: center; padding: 12px 16px; }
    .submit { grid-row: 4; grid-column: 2; display: flex; align-items: center; justify-content: center; }

    .map {
      width: 100%;
      height: calc(100% - 42px);
      border-radius: 0;
      border: 1.8px solid #111827;
      background:
        linear-gradient(rgba(30,144,255,0.22) 1px, transparent 1px),
        linear-gradient(90deg, rgba(30,144,255,0.22) 1px, transparent 1px),
        #f8fafc;
      background-size: 36px 36px;
      position: relative;
    }
    .region { position: absolute; border-radius: 0; opacity: 0.82; border: 2px solid rgba(0,0,0,0.55); }
    .r-blue { background: var(--blue); }
    .r-green { background: var(--green); }
    .r-yellow { background: var(--yellow); }
    .r-purple { background: var(--purple); }
    .r-brown { background: var(--brown); }
    .obs { background: #111827; opacity: 0.72; border-radius: 0; }
    .trajectory-layer {
      position: absolute;
      inset: 0;
      pointer-events: none;
      z-index: 8;
    }
    .trajectory-layer svg {
      width: 100%;
      height: 100%;
      overflow: visible;
    }
    .trajectory-layer polyline {
      fill: none;
      stroke: #0ea5e9;
      stroke-width: 0.75;
      stroke-linejoin: round;
      stroke-linecap: round;
    }
    .vehicle {
      position: absolute;
      width: 18px; height: 18px;
      transform: translate(-50%, -50%);
      border-radius: 50%;
      background: #dc2626;
      border: 3px solid white;
      box-shadow: 0 0 0 1px #991b1b;
      z-index: 10;
    }

    .formula-box {
      min-height: 76px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      font-size: 16px;
      border: 1px solid #d1d5db;
      border-radius: 10px;
      background: #f9fafb;
      padding: 12px;
      text-align: center;
    }
    .automaton-canvas {
      height: calc(100% - 32px);
      border-radius: 10px;
      border: 1px solid #e5e7eb;
      display: flex;
      align-items: flex-start;
      justify-content: center;
      color: var(--muted);
      background: #fbfdff;
      padding: 0 2px 2px 2px;
      overflow: hidden;
    }
    .quotient-svg {
      width: 100%;
      height: 118%;
      transform: translateY(-9%);
    }
    .legend {
      display: grid;
      grid-template-columns: 24px 1fr 24px 1fr;
      gap: 8px 10px;
      align-items: center;
      font-size: 14px;
    }
    .swatch { width: 24px; height: 18px; border: 1.5px solid #111827; border-radius: 2px; }
    .info-text {
      white-space: pre-wrap;
      font-size: 15px;
      line-height: 1.45;
      color: #1f2937;
      overflow-y: auto;
      max-height: calc(100% - 84px);
      padding-right: 6px;
    }
    .validation-controls {
      display: none;
      gap: 10px;
      margin-top: 16px;
    }
    .validation-controls.active {
      display: grid;
      grid-template-columns: 1fr 1fr;
    }
    .validation-controls button {
      height: 46px;
      font-size: 16px;
      border-radius: 10px;
    }
    #validate-yes {
      background: #15803d;
    }
    #validate-no {
      background: #991b1b;
    }
    #validate-yes:hover {
      background: #166534;
    }
    #validate-no:hover {
      background: #7f1d1d;
    }
    textarea {
      width: 100%; height: 100%; resize: none;
      border: none; outline: none;
      font-size: 18px;
      font-family: inherit;
    }
    button {
      width: 100%; height: 100%;
      border: none;
      border-radius: 12px;
      background: #111827;
      color: white;
      font-size: 22px;
      font-weight: 800;
      cursor: pointer;
    }
    button:hover { background: #374151; }
  </style>
</head>
<body>
  <main class="app">
    <section class="card workspace">
      <div class="title">Abstract system</div>
      <div class="map">
        $(genie_workspace_map_html())
      </div>
    </section>

    <section class="card formula">
      <div class="title">LTL Formula <span style="font-weight:400; font-size:16px; color:#4b5563;">(G(!obstacle) & ...)</span></div>
      <div id="formula" class="formula-box">Formula: &lt;none&gt;</div>
    </section>



    <section class="card ap">
      <div class="title">Atomic propositions</div>
      <div class="legend">
        <div class="swatch" style="background: var(--blue)"></div><div>blue</div>
        <div class="swatch" style="background: var(--green)"></div><div>green</div>
        <div class="swatch" style="background: var(--yellow)"></div><div>yellow</div>
        <div class="swatch" style="background: var(--purple)"></div><div>purple</div>
        <div class="swatch" style="background: var(--brown)"></div><div>brown</div>
        <div class="swatch" style="background:#111827"></div><div>obstacle</div>
      </div>
    </section>

    <section class="card info">
      <div class="title">Info</div>
      <div id="info" class="info-text">Status: waiting for requirement</div>
      <div id="validation-controls" class="validation-controls">
        <button id="validate-yes" type="button">Yes, correct</button>
        <button id="validate-no" type="button">No, rephrase</button>
      </div>
    </section>

    <section class="card input">
      <textarea id="requirement" placeholder="You can enter your requirements here..."></textarea>
    </section>

    <section class="card submit">
      <button id="submit">Submit</button>
    </section>
  </main>

  <script>
    const button = document.getElementById('submit');
    const req = document.getElementById('requirement');
    const formula = document.getElementById('formula');
    const info = document.getElementById('info');
    const automaton = document.getElementById('automaton');
    const vehicle = document.getElementById('vehicle');
    const trajectoryPolyline = document.getElementById('trajectory-polyline');
    const validationControls = document.getElementById('validation-controls');
    const validateYes = document.getElementById('validate-yes');
    const validateNo = document.getElementById('validate-no');
    const quotientSvgHtml = automaton ? automaton.innerHTML : '';
    let animationToken = 0;
    let pendingPlanToken = null;
    let currentWorkspacePoint = null;

    function cancelCurrentAnimation() {
      animationToken += 1;
      trajectoryPolyline.setAttribute('points', '');
      setActiveQuotientState('q0');
    }

    function workspaceToSvgPoint(point) {
      const x = 10 * point[0];
      const y = 100 - 10 * point[1];
      return x.toFixed(2) + ',' + y.toFixed(2);
    }

    function workspaceToCssPoint(point) {
      const left = 10 * point[0];
      const top = 100 - 10 * point[1];
      return { left: left.toFixed(2) + '%', top: top.toFixed(2) + '%' };
    }

    function setVehicleWorkspacePoint(point) {
      if (!Array.isArray(point) || point.length < 2) return;
      currentWorkspacePoint = [Number(point[0]), Number(point[1])];
      const css = workspaceToCssPoint(currentWorkspacePoint);
      vehicle.style.left = css.left;
      vehicle.style.top = css.top;
    }

    function sleep(ms) {
      return new Promise(resolve => setTimeout(resolve, ms));
    }

    function setValidationControlsActive(active) {
      validationControls.classList.toggle('active', Boolean(active));
      validateYes.disabled = !active;
      validateNo.disabled = !active;
    }

    async function postJson(url, payload) {
      const response = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload || {})
      });
      const data = await response.json();
      if (!data.ok) throw new Error(data.error || 'Unknown error');
      return data;
    }

    function showValidationDraft(data) {
      pendingPlanToken = (data.plan_token === undefined || data.plan_token === null) ? null : data.plan_token;
      formula.textContent = data.formula || 'Formula: <none>';
      if (automaton) {
        automaton.innerHTML = quotientSvgHtml;
        setActiveQuotientState('q0');
      }
      trajectoryPolyline.setAttribute('points', '');
      let text = data.validation_message || data.status || 'Please validate the generated requirement.';
      text = text.replace(/^\\s*\\*\\*?Validation Message:?\\*\\*?\\s*/i, '');
      text = text.replace(/^\\s*Validation Message:?\\s*/i, '');
      if (data.backtranslation) {
        text += String.fromCharCode(10) + String.fromCharCode(10) + 'Backtranslation: ' + data.backtranslation;
      }
      if (data.quotient_message && String(data.quotient_message).trim() !== 'Requirement passed the coarse AP-label quotient compatibility check.') {
        text += String.fromCharCode(10) + String.fromCharCode(10) + data.quotient_message;
      }
      info.textContent = text;
      setValidationControlsActive(true);
    }

    function clearActiveQuotientState() {
      if (!automaton) return;
      document.querySelectorAll('.qnode').forEach(node => node.classList.remove('qactive'));
    }

    function setActiveQuotientState(qid) {
      if (!automaton) return;
      clearActiveQuotientState();
      const node = document.getElementById('quotient-' + qid);
      if (node) node.classList.add('qactive');
    }

    async function animateTrajectory(points, dtSeconds, quotientStates, localAnimationToken) {
      trajectoryPolyline.setAttribute('points', '');
      if (!Array.isArray(points) || points.length === 0) return;
      if (localAnimationToken !== animationToken) return;

      const stepMs = Math.max(1, Math.round(1000 * (Number(dtSeconds) || 0.3)));
      const drawn = [];

      if (points.length === 1) {
        setVehicleWorkspacePoint(points[0]);
        if (Array.isArray(quotientStates) && quotientStates.length > 0) {
          setActiveQuotientState(quotientStates[0]);
        }
        return;
      }

      for (let i = 0; i < points.length; i++) {
        if (localAnimationToken !== animationToken) return;
        const point = points[i];
        drawn.push(point);
        trajectoryPolyline.setAttribute('points', drawn.map(workspaceToSvgPoint).join(' '));
        setVehicleWorkspacePoint(point);
        if (Array.isArray(quotientStates) && quotientStates[i]) {
          setActiveQuotientState(quotientStates[i]);
        }
        await sleep(stepMs);
        if (localAnimationToken !== animationToken) return;
      }
    }

    button.addEventListener('click', async () => {
      console.log('Submit clicked');
      const nl = req.value.trim();
      if (!nl) {
        info.textContent = 'Status: please enter a requirement first.';
        return;
      }

      cancelCurrentAnimation();
      pendingPlanToken = null;
      info.textContent = 'Previous execution interrupted. Preparing the new requirement...';
      button.disabled = true;
      button.textContent = 'Checking...';
      setValidationControlsActive(false);
      info.textContent = 'Status: translating requirement and preparing validation question...';
      formula.textContent = 'Formula: ...';
      if (automaton) {
        automaton.innerHTML = '<div style="padding: 18px; color: #6b7280;">Preparing validation...</div>';
      }

      try {
        const data = await postJson('/api/submit', {
          requirement: nl,
          current_state: currentWorkspacePoint
        });
        showValidationDraft(data);
      } catch (err) {
        formula.textContent = 'Formula: <error>';
        if (automaton) {
          automaton.innerHTML = quotientSvgHtml;
          setActiveQuotientState('q0');
        }
        if (err instanceof TypeError && String(err.message).toLowerCase().includes('fetch')) {
          info.textContent = 'Failed to reach the Genie server. The Julia process probably stopped or crashed. Check the Julia terminal for the first error above the shutdown message.';
        } else {
          info.textContent = err.message;
        }
        trajectoryPolyline.setAttribute('points', '');
        console.error(err);
      } finally {
        button.disabled = false;
        button.textContent = 'Submit';
      }
    });

    validateYes.addEventListener('click', async () => {
      const acceptedPlanToken = pendingPlanToken;
      cancelCurrentAnimation();
      setValidationControlsActive(false);
      validateYes.disabled = true;
      validateNo.disabled = true;
      info.textContent = 'Validation accepted. Planning...';
      button.disabled = true;
      button.textContent = 'Planning...';
      try {
        const data = await postJson('/api/confirm', {});
        if (acceptedPlanToken !== null && data.plan_token !== undefined && data.plan_token !== acceptedPlanToken) {
          info.textContent = 'A newer requirement replaced this one. Ignoring the stale plan.';
          return;
        }
        formula.textContent = data.formula || 'Formula: <none>';
        const ntraj = (data.trajectory || []).length;
        if (automaton) {
          automaton.innerHTML = quotientSvgHtml;
          setActiveQuotientState('q0');
        }
        if (ntraj <= 1) {
          info.textContent = (data.status || 'Done.') + String.fromCharCode(10) + 'No nontrivial trajectory was returned by the runtime.';
        } else {
          info.textContent = data.status || 'Done.';
        }
        button.disabled = false;
        button.textContent = 'Submit';
        const localAnimationToken = animationToken;
        animateTrajectory(data.trajectory || [], data.dt_seconds || 0.3, data.quotient_states || [], localAnimationToken)
          .catch(err => {
            if (localAnimationToken === animationToken) {
              info.textContent = err.message;
              console.error(err);
            }
          });
      } catch (err) {
        formula.textContent = 'Formula: <error>';
        if (automaton) {
          automaton.innerHTML = quotientSvgHtml;
          setActiveQuotientState('q0');
        }
        info.textContent = err.message;
        trajectoryPolyline.setAttribute('points', '');
        console.error(err);
      } finally {
        if (button.textContent === 'Planning...') {
          button.disabled = false;
          button.textContent = 'Submit';
        }
      }
    });

    validateNo.addEventListener('click', async () => {
      cancelCurrentAnimation();
      pendingPlanToken = null;
      setValidationControlsActive(false);
      try {
        const data = await postJson('/api/reject', {});
        info.textContent = data.status || 'Please rephrase the requirement and submit again.';
        formula.textContent = 'Formula: <none>';
        trajectoryPolyline.setAttribute('points', '');
        if (automaton) {
          automaton.innerHTML = quotientSvgHtml;
          setActiveQuotientState('q0');
        }
        req.focus();
      } catch (err) {
        info.textContent = err.message;
        console.error(err);
      }
    });
  </script>
</body>
</html>
"""
end

function setup_genie_routes!()
    Genie.Router.route("/", method = GET) do
        return Genie.Renderer.Html.html(genie_index_html())
    end

    Genie.Router.route("/api/submit", method = POST) do
        try
            payload_text = String(rawpayload())
            isempty(strip(payload_text)) && error("Empty request body.")

            payload = JSON3.read(payload_text)
            nl = haskey(payload, :requirement) ? String(payload.requirement) : ""
            current_state_xy = haskey(payload, :current_state) ? payload.current_state : nothing
            result = genie_prepare_validation(nl, current_state_xy)
            return Genie.Renderer.Json.json(result)
        catch err
            msg = sprint(showerror, err)
            bt = catch_backtrace()
            @error "Genie /api/submit failed" exception=(err, bt)

            lower_msg = lowercase(msg)
            user_msg = msg
            if occursin("coarse ap-label quotient", lower_msg) ||
               occursin("no reachable accepting", lower_msg) ||
               occursin("quotient compatibility check failed", lower_msg) ||
               occursin("not compatible", lower_msg)
                user_msg = "Requirement is not compatible with the current quotient abstraction: no reachable accepting Büchi run exists in the quotient product."
            elseif occursin("quotient", lower_msg) ||
                   occursin("infeasible", lower_msg) ||
                   occursin("no_admissible_control", lower_msg)
                user_msg = "Request failed during quotient precheck or controller synthesis. Please try a weaker or differently phrased requirement."
            elseif occursin("api key", lower_msg) || occursin("authentication", lower_msg) || occursin("401", lower_msg)
                user_msg = "The language-model API call failed because of a missing or invalid API key. Check ANTHROPIC_API_KEY / OPENAI_API_KEY in the Julia environment before starting Genie."
            elseif occursin("connection", lower_msg) || occursin("http", lower_msg) || occursin("request", lower_msg)
                user_msg = "A planner or language-model HTTP request failed. Check the Julia terminal for the detailed error log."
            end

            return Genie.Renderer.Json.json(Dict(
                "ok" => false,
                "error" => user_msg,
                "debug_error" => msg,
            ))
        end
    end

    Genie.Router.route("/api/confirm", method = POST) do
        try
            result = genie_confirm_validation_and_plan()
            return Genie.Renderer.Json.json(result)
        catch err
            msg = sprint(showerror, err)
            bt = catch_backtrace()
            @error "Genie /api/confirm failed" exception=(err, bt)
            lower_msg = lowercase(msg)
            user_msg = msg
            if occursin("coarse ap-label quotient", lower_msg) ||
               occursin("no reachable accepting", lower_msg) ||
               occursin("quotient compatibility check failed", lower_msg) ||
               occursin("not compatible", lower_msg)
                user_msg = "Requirement is not compatible with the current quotient abstraction: no reachable accepting Büchi run exists in the quotient product."
            elseif occursin("quotient", lower_msg) ||
                   occursin("infeasible", lower_msg) ||
                   occursin("no_admissible_control", lower_msg)
                user_msg = "Request failed during quotient precheck or controller synthesis. Please try a weaker or differently phrased requirement."
            end
            return Genie.Renderer.Json.json(Dict(
                "ok" => false,
                "error" => user_msg,
                "debug_error" => msg,
            ))
        end
    end

    Genie.Router.route("/api/reject", method = POST) do
        try
            result = genie_reject_validation()
            return Genie.Renderer.Json.json(result)
        catch err
            msg = sprint(showerror, err)
            bt = catch_backtrace()
            @error "Genie /api/reject failed" exception=(err, bt)
            return Genie.Renderer.Json.json(Dict(
                "ok" => false,
                "error" => msg,
            ))
        end
    end

    return nothing
end

function start_genie_interface(; host::AbstractString = GENIE_DEFAULT_HOST,
                                port::Integer = GENIE_DEFAULT_PORT,
                                async::Bool = true)
    Genie.config.run_as_server = true
    Genie.config.server_host = host
    Genie.config.server_port = port
    setup_genie_routes!()
    println("Genie interface available at http://$(host):$(port)")
    return Genie.up(port, host; async = async)
end
