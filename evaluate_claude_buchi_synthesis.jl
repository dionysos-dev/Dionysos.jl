#!/usr/bin/env julia

"""
Evaluate controller synthesis for the correctly translated formulas in
`Claude_task_specific_rewritten.json`.

The script:
  1. keeps records with `equivalent == true` and a nonempty `predicted_ltl`;
  2. translates each unique formula to the deterministic, state-based Büchi
     representation required by the existing solver;
  3. calls `synthesize_buchi_controller(...; semantics=:robust)`;
  4. consumes each initial plant label before checking initial product states;
  5. verifies universal infinite-horizon Büchi acceptance of the closed-loop
     product under all nondeterministic plant successors.

The verifier explores through accepting states.  Satisfaction requires a
nonblocking closed loop with no undefined transition and no reachable cycle made
entirely of non-accepting Büchi states.

Usage:
    julia --project=. evaluate_claude_buchi_synthesis.jl
    julia --project=. evaluate_claude_buchi_synthesis.jl INPUT.json OUTPUT.json
    julia --project=. evaluate_claude_buchi_synthesis.jl INPUT.json OUTPUT.json 100

The optional third argument limits the number of eligible dataset records, which
is useful for a smoke test.  Exact duplicate formulas are synthesized only once,
but the final rates are record-weighted and every eligible record is retained.
"""

using Dates
using JSON3
using CSV
using DataFrames
using JLD2
using StaticArrays
using Dionysos

const EVAL_ROOT = @__DIR__
const DEFAULT_INPUT = joinpath(EVAL_ROOT, "Claude_task_specific_rewritten.json")
const DEFAULT_OUTPUT = joinpath(EVAL_ROOT, "Claude_task_specific_buchi_synthesis_results.json")
const DEFAULT_ABSTRACTION = joinpath(EVAL_ROOT, "abstraction.jld2")
const TASK_APS_EVAL = Set(["obs", "blue", "green", "purple", "brown", "yellow"])
const EVALUATION_SEMANTICS_VERSION = "deterministic-buchi-v1-initial-label-infinite-verification"

if !isdefined(Main, :LTL2Automata)
    include(joinpath(EVAL_ROOT, "src", "LTL", "LTL2Automata.jl"))
end
const EVAL_L2A = Main.LTL2Automata

field_or(record, name::Symbol, default=nothing) = haskey(record, name) ? getproperty(record, name) : default

function canonical_formula(formula::String)
    try
        return string(EVAL_L2A.Spot.SpotFormula(formula))
    catch
        return strip(formula)
    end
end

string_dict(object) = Dict{String,Any}(String(key) => value for (key, value) in pairs(object))

"""Load only the saved symbolic abstraction; no Genie/GUI/LLM code is included."""
function load_saved_abstract_system(path::String)
    isfile(path) || error(
        "Saved abstraction not found: $path. Set BUCHI_ABSTRACTION_PATH or create abstraction.jld2 first."
    )
    println("Loading symbolic abstraction from: ", path)
    return JLD2.jldopen(path, "r"; mmaparrays=true) do file
        haskey(file, "abstract_system") || error("$path has no `abstract_system` entry")
        read(file, "abstract_system")
    end
end

"""
Construct the exact TS and AP labels needed by the Büchi solver.

This is the non-interface subset of the setup in `InteractiveWithFeedBack.jl`.
"""
function build_batch_backend(abstraction_path::String=DEFAULT_ABSTRACTION)
    abstract_system = load_saved_abstract_system(abstraction_path)
    UT = Dionysos.Utils
    SY = Dionysos.Symbolic
    DO = Dionysos.Domain

    theta_min = -pi - 0.4
    theta_max = pi + 0.4
    rect(x1, y1, x2, y2) = UT.HyperRectangle(
        SVector(x1, y1, theta_min), SVector(x2, y2, theta_max)
    )

    obstacle_margin = 0.12
    x1_lb = [0.5, 0.5, 3.5, 3.5, 4.0, 5.5]
    x1_ub = [2.5, 2.5, 4.0, 4.0, 6.0, 6.0]
    x2_lb = [5.0, 6.5, 1.0, 4.0, 4.0, 4.5]
    x2_ub = [5.5, 9.0, 3.0, 6.0, 4.5, 6.0]
    obstacles = [
        rect(
            max(0.0, x1_lb[i] - obstacle_margin),
            max(0.0, x2_lb[i] - obstacle_margin),
            min(10.0, x1_ub[i] + obstacle_margin),
            min(10.0, x2_ub[i] + obstacle_margin),
        ) for i in eachindex(x1_lb)
    ]
    regions = Dict{Symbol,Any}(
        :obs => obstacles,
        :blue => rect(6.0, 7.5, 9.0, 9.5),
        :green => rect(5.0, 1.0, 7.5, 2.5),
        :purple => rect(7.0, 3.5, 8.5, 6.5),
        :brown => rect(4.0, 7.0, 7.5, 8.0),
        :yellow => rect(0.5, 1.0, 2.5, 3.5),
    )
    labeler, _ = EVAL_L2A.labeler_and_inner_sets(
        abstract_system; regions=regions, forbidden_aps=[:obs]
    )

    x0_cont = SVector(4.8, 5.0, pi / 2)
    hx = SVector(0.20, 0.20, 0.30)
    start_region = UT.HyperRectangle(
        x0_cont - SVector(0.10, 0.10, 0.10),
        x0_cont + SVector(0.10, 0.10, 0.10),
    )
    raw_ts = EVAL_L2A.ts_from_symbolic_onepass(
        abstract_system, start_region; reachable_only=false, x0_abs=nothing
    )

    initial_cell = UT.HyperRectangle(x0_cont - hx / 2, x0_cont + hx / 2)
    initial_states = SY.get_states_from_set(abstract_system, initial_cell, DO.CENTER)
    isempty(initial_states) &&
        (initial_states = SY.get_states_from_set(abstract_system, initial_cell, DO.OUTER))
    isempty(initial_states) && error("Could not map the configured initial state into the abstraction")

    ts = EVAL_L2A.TSFromAdj(raw_ts.Adj, [first(initial_states)])
    labels_ts, _ = EVAL_L2A.precompute_labels_and_forbidden(ts, labeler, Symbol[])
    println("Batch backend ready: $(EVAL_L2A.nstates(ts)) states, $(EVAL_L2A.ninputs(ts)) inputs")
    return (; abstract_system, ts, labeler, labels_ts, x0_abs_start=first(initial_states))
end

q_to_index(aut, q::Int) = aut.zero_based ? q + 1 : q
index_to_q(aut, iq::Int) = aut.zero_based ? iq - 1 : iq

product_id(nx::Int, x::Int, iq::Int) = (iq - 1) * nx + x

function product_components(nx::Int, pid::Int)
    iq = div(pid - 1, nx) + 1
    x = mod(pid - 1, nx) + 1
    return x, iq
end

"""Return true iff a finite directed graph contains a cycle."""
function has_directed_cycle(nodes::Vector{Int}, adjacency::Dict{Int,Vector{Int}})
    isempty(nodes) && return false
    node_set = Set(nodes)
    indegree = Dict(node => 0 for node in nodes)

    for src in nodes
        for dst in get(adjacency, src, Int[])
            dst in node_set || continue
            indegree[dst] = indegree[dst] + 1
        end
    end

    queue = [node for node in nodes if indegree[node] == 0]
    head = 1
    removed = 0
    while head <= length(queue)
        src = queue[head]
        head += 1
        removed += 1
        for dst in get(adjacency, src, Int[])
            dst in node_set || continue
            indegree[dst] -= 1
            indegree[dst] == 0 && push!(queue, dst)
        end
    end
    return removed != length(nodes)
end

"""
Verify an infinite-horizon Büchi policy against all plant successors.

The initial plant label is consumed before entering the product.  A reachable
undefined automaton transition, missing control, plant dead end, or cycle wholly
outside the accepting set makes the check fail.
"""
function verify_closed_loop_buchi(ts, aut, policy, next_q, initial_states)
    nx = EVAL_L2A.nstates(ts)
    iq0 = q_to_index(aut, aut.init)
    starts = Int[]
    violation_count = 0
    for x0 in initial_states
        iq_after = next_q[iq0, x0]
        if iq_after == EVAL_L2A.I_INVALID_Q
            violation_count += 1
        else
            push!(starts, product_id(nx, x0, iq_after))
        end
    end

    queue = copy(starts)
    visited = Set(starts)
    adjacency = Dict{Int,Vector{Int}}()
    head = 1

    while head <= length(queue)
        pid = queue[head]
        head += 1
        x, iq = product_components(nx, pid)

        if pid < 1 || pid > length(policy)
            violation_count += 1
            continue
        end

        u = policy[pid]
        if u == 0
            violation_count += 1
            continue
        end

        successors = EVAL_L2A.post(ts, x, u)
        if isempty(successors)
            violation_count += 1
            continue
        end

        product_successors = get!(adjacency, pid, Int[])
        for xp in successors
            iqp = next_q[iq, xp]
            if iqp == EVAL_L2A.I_INVALID_Q
                violation_count += 1
                continue
            end
            pidp = product_id(nx, xp, iqp)
            push!(product_successors, pidp)
            if !(pidp in visited)
                push!(visited, pidp)
                push!(queue, pidp)
            end
        end
        unique!(product_successors)
    end

    nonaccepting = Int[]
    for pid in visited
        _, iq = product_components(nx, pid)
        q = index_to_q(aut, iq)
        EVAL_L2A.is_accepting(aut, q) || push!(nonaccepting, pid)
    end
    nonaccepting_cycle = has_directed_cycle(nonaccepting, adjacency)
    satisfied = violation_count == 0 && !nonaccepting_cycle

    return Dict(
        "satisfied" => satisfied,
        "reachable_product_states" => length(visited),
        "nonaccepting_product_states" => length(nonaccepting),
        "has_reachable_nonaccepting_cycle" => nonaccepting_cycle,
        "transition_violation_count" => violation_count,
    )
end

function atomic_write_json(path::String, value)
    mkpath(dirname(path))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
    mv(temporary, path; force=true)
    return path
end

function automaton_signature(aut)
    transitions = sort!(collect(aut.delta); by=pair -> (pair.first[1], pair.first[2], pair.second))
    transition_tuple = Tuple((pair.first[1], pair.first[2], pair.second) for pair in transitions)
    return (Tuple(aut.aps), aut.init, aut.nstates, Tuple(aut.accepting), transition_tuple)
end

function rate(numerator::Integer, denominator::Integer)
    denominator == 0 && return nothing
    return numerator / denominator
end

function evaluate_formula(
    formula::String,
    backend,
    labels_ts;
    max_steps::Int=1200,
    automaton_result_cache::Union{Nothing,AbstractDict}=nothing,
)
    result = Dict{String,Any}(
        "formula" => formula,
        "synthesis_completed" => false,
        "controller_realizable" => false,
        "initial_state_feasible" => false,
        "specification_satisfied" => false,
        "error" => nothing,
    )

    try
        aut = nothing
        translation_seconds = @elapsed aut = EVAL_L2A.translate_ltl_buchi(formula)

        unknown_aps = sort!(collect(setdiff(Set(aut.aps), TASK_APS_EVAL)))
        isempty(unknown_aps) || error("Formula contains APs without environment labels: $(join(unknown_aps, ", "))")

        signature = automaton_signature(aut)
        if automaton_result_cache !== nothing && haskey(automaton_result_cache, signature)
            reused = copy(automaton_result_cache[signature])
            reused["formula"] = formula
            reused["translation_seconds"] = translation_seconds
            reused["reused_automaton_result"] = true
            return reused
        end

        WinP = polP = next_qP = infoP = nothing
        synthesis_seconds = @elapsed begin
            WinP, polP, next_qP, infoP = EVAL_L2A.synthesize_buchi_controller(
                backend.ts, aut, labels_ts; semantics=:robust
            )
        end

        nx = EVAL_L2A.nstates(backend.ts)
        iq0 = q_to_index(aut, aut.init)
        initial_product_ids = Int[]
        initial_transition_valid = true
        for x0 in EVAL_L2A.initset(backend.ts)
            iq_after = next_qP[iq0, x0]
            if iq_after == EVAL_L2A.I_INVALID_Q
                initial_transition_valid = false
            else
                push!(initial_product_ids, product_id(nx, x0, iq_after))
            end
        end
        controller_realizable = any(WinP)
        initial_state_feasible =
            initial_transition_valid && !isempty(initial_product_ids) &&
            all(pid -> WinP[pid], initial_product_ids)

        # Record synthesis outcomes immediately.  If the later simulation or
        # satisfaction check fails, synthesis itself still completed.
        merge!(result, Dict(
            "synthesis_completed" => true,
            "controller_realizable" => controller_realizable,
            "initial_state_feasible" => initial_state_feasible,
            "translation_seconds" => translation_seconds,
            "synthesis_seconds" => synthesis_seconds,
            "automaton_states" => aut.nstates,
            "automaton_accepting_states" => length(aut.accepting),
            "winning_product_states" => count(WinP),
            "total_product_states" => length(WinP),
        ))

        simulation_status = "not_run"
        simulation_steps = 0
        certificate = Dict{String,Any}(
            "satisfied" => false,
            "reachable_product_states" => 0,
            "nonaccepting_product_states" => 0,
            "has_reachable_nonaccepting_cycle" => false,
            "transition_violation_count" => 0,
        )

        if initial_state_feasible
            certificate = verify_closed_loop_buchi(
                backend.ts, aut, polP, next_qP, EVAL_L2A.initset(backend.ts)
            )
            xs, _, _, status = EVAL_L2A.simulate_buchi_policy(
                backend.ts, aut, labels_ts, polP, next_qP;
                x0=first(EVAL_L2A.initset(backend.ts)), q0=aut.init,
                max_steps=max_steps, stop_at_accept=false,
                consume_initial_label=true,
            )
            simulation_status =
                certificate["satisfied"] && status == :max_steps ? "persistent" : string(status)
            simulation_steps = max(0, length(xs) - 1)
        end

        merge!(result, Dict(
            "specification_satisfied" => initial_state_feasible && certificate["satisfied"],
            "representative_simulation_status" => simulation_status,
            "representative_simulation_steps" => simulation_steps,
            "verification_reachable_product_states" => certificate["reachable_product_states"],
            "verification_nonaccepting_product_states" => certificate["nonaccepting_product_states"],
            "verification_has_nonaccepting_cycle" => certificate["has_reachable_nonaccepting_cycle"],
            "verification_transition_violation_count" => certificate["transition_violation_count"],
            "reused_automaton_result" => false,
        ))
        automaton_result_cache !== nothing &&
            (automaton_result_cache[signature] = copy(result))
    catch err
        result["error"] = sprint(showerror, err)
    end
    return result
end

function make_summary(outcomes)
    attempted = length(outcomes)
    completed = count(r -> r["synthesis_completed"] === true, outcomes)
    realizable = count(r -> r["controller_realizable"] === true, outcomes)
    initial_feasible = count(r -> r["initial_state_feasible"] === true, outcomes)
    satisfied = count(r -> r["specification_satisfied"] === true, outcomes)
    return Dict(
        "eligible_records" => attempted,
        "synthesis_completed_records" => completed,
        "controller_realizable_records" => realizable,
        "initial_state_feasible_records" => initial_feasible,
        "specification_satisfied_records" => satisfied,
        "synthesis_completion_rate" => rate(completed, attempted),
        "controller_realizability_rate_among_completed" => rate(realizable, completed),
        "initial_state_feasibility_rate_among_completed" => rate(initial_feasible, completed),
        "initial_state_feasibility_rate_among_realizable" => rate(initial_feasible, realizable),
        "specification_satisfaction_rate_among_initial_state_feasible" => rate(satisfied, initial_feasible),
        "end_to_end_satisfied_rate" => rate(satisfied, attempted),
    )
end

function write_csv(path::String, records)
    rows = [(
        record_id = r["record_id"],
        input_field = r["input_field"],
        source_text = r["source_text"],
        original_ltl = r["original_ltl"],
        predicted_ltl = r["predicted_ltl"],
        synthesis_completed = r["synthesis_completed"],
        controller_realizable = r["controller_realizable"],
        initial_state_feasible = r["initial_state_feasible"],
        specification_satisfied = r["specification_satisfied"],
        automaton_states = get(r, "automaton_states", missing),
        winning_product_states = get(r, "winning_product_states", missing),
        total_product_states = get(r, "total_product_states", missing),
        translation_seconds = get(r, "translation_seconds", missing),
        synthesis_seconds = get(r, "synthesis_seconds", missing),
        representative_simulation_status = get(r, "representative_simulation_status", missing),
        error = something(get(r, "error", nothing), ""),
    ) for r in records]
    CSV.write(path, DataFrame(rows))
    return path
end

function main(args=ARGS)
    input_path = length(args) >= 1 ? abspath(args[1]) : DEFAULT_INPUT
    output_path = length(args) >= 2 ? abspath(args[2]) : DEFAULT_OUTPUT
    limit = length(args) >= 3 ? parse(Int, args[3]) : typemax(Int)
    isfile(input_path) || error("Input dataset not found: $input_path")

    dataset = JSON3.read(read(input_path, String))
    eligible = [record for record in dataset.results if
        field_or(record, :equivalent, false) === true &&
        field_or(record, :predicted_ltl, nothing) !== nothing &&
        !isempty(strip(String(record.predicted_ltl)))
    ]
    eligible = eligible[1:min(length(eligible), limit)]
    exact_formulas = unique(String(record.predicted_ltl) for record in eligible)
    key_to_formula = Dict{String,String}()
    formula_keys = String[]
    for formula in exact_formulas
        key = canonical_formula(formula)
        if !haskey(key_to_formula, key)
            key_to_formula[key] = formula
            push!(formula_keys, key)
        end
    end

    println("Eligible records: ", length(eligible))
    println("Exact unique formulas:     ", length(exact_formulas))
    println("Canonical formula classes: ", length(formula_keys))
    abstraction_path = abspath(get(ENV, "BUCHI_ABSTRACTION_PATH", DEFAULT_ABSTRACTION))
    backend = build_batch_backend(abstraction_path)
    labels_ts = backend.labels_ts

    cache = Dict{String,Dict{String,Any}}()
    automaton_result_cache = Dict{Any,Dict{String,Any}}()
    checkpoint_path = output_path * ".checkpoint.json"
    if isfile(checkpoint_path)
        checkpoint = JSON3.read(read(checkpoint_path, String))
        if field_or(checkpoint, :semantics_version, "") == EVALUATION_SEMANTICS_VERSION &&
           field_or(checkpoint, :input_path, "") == input_path
            for prior in checkpoint.formula_results
                prior_dict = string_dict(prior)
                key = String(prior_dict["canonical_formula"])
                key in formula_keys && (cache[key] = prior_dict)
            end
            println("Resumed ", length(cache), " canonical formula results from checkpoint.")
        end
    end

    processed = 0
    for (index, key) in enumerate(formula_keys)
        if haskey(cache, key)
            processed += 1
            continue
        end
        formula = key_to_formula[key]
        println("[$index/$(length(formula_keys))] ", formula)
        outcome = evaluate_formula(
            formula, backend, labels_ts; automaton_result_cache=automaton_result_cache
        )
        outcome["canonical_formula"] = key
        cache[key] = outcome
        processed += 1

        if processed % 25 == 0
            checkpoint = Dict(
                "status" => "in_progress",
                "semantics_version" => EVALUATION_SEMANTICS_VERSION,
                "input_path" => input_path,
                "completed_canonical_formulas" => length(cache),
                "total_canonical_formulas" => length(formula_keys),
                "formula_results" => collect(values(cache)),
            )
            atomic_write_json(checkpoint_path, checkpoint)
        end
    end

    records = Dict{String,Any}[]
    for record in eligible
        formula = String(record.predicted_ltl)
        key = canonical_formula(formula)
        outcome = copy(cache[key])
        merge!(outcome, Dict(
            "record_id" => field_or(record, :record_id),
            "input_field" => String(field_or(record, :input_field, "")),
            "source_model" => String(field_or(record, :source_model, "")),
            "source_text" => String(field_or(record, :source_text, "")),
            "original_ltl" => String(field_or(record, :original_ltl, "")),
            "predicted_ltl" => formula,
            "formula" => formula,
            "canonical_formula" => key,
            "equivalent" => true,
        ))
        push!(records, outcome)
    end

    unique_outcomes = collect(values(cache))
    report = Dict(
        "created_at" => string(now()),
        "semantics_version" => EVALUATION_SEMANTICS_VERSION,
        "input_path" => input_path,
        "methodology" => Dict(
            "filter" => "equivalent == true and predicted_ltl is nonempty",
            "solver" => "synthesize_buchi_controller",
            "game_semantics" => "robust: exists control, forall plant successors",
            "synthesis_completed" => "the solver returned WinP and a policy without throwing an error",
            "controller_realizable" => "WinP is nonempty: at least one product state is winning",
            "initial_state_feasible" => "all configured initial product states belong to WinP",
            "specification_satisfaction" => "the closed-loop product is nonblocking and every infinite path visits accepting states infinitely often",
            "important_note" => "Initial plant labels are consumed before product membership is checked; representative simulation is illustrative only.",
        ),
        "record_weighted_summary" => make_summary(records),
        "unique_formula_summary" => merge(
            make_summary(unique_outcomes),
            Dict(
                "exact_unique_formulas" => length(exact_formulas),
                "canonical_formula_classes" => length(unique_outcomes),
            ),
        ),
        "records" => records,
    )

    atomic_write_json(output_path, report)
    csv_path = splitext(output_path)[1] * ".csv"
    write_csv(csv_path, records)
    println("Wrote JSON report: ", output_path)
    println("Wrote CSV table:  ", csv_path)
    println(JSON3.pretty(report["record_weighted_summary"]))
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
