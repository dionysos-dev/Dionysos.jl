using StaticArrays, JuMP, Plots, MathematicalSystems, HybridSystems, GLMakie, Serialization
using LinearAlgebra: norm
using Dionysos
import MathOptInterface as MOI



const DI = Dionysos
const UT = DI.Utils
const DO = DI.Utils
const ST = DI.System
const SY = DI.Symbolic
const PR = DI.Problem
const OP = DI.Optim
const AB = OP.Abstraction
const PB = Dionysos.Problem

# Helper to construct a GridFree object from a grid spacing vector.
# In the current Dionysos API, GridFree is defined in the mapping module
# used internally by UniformGridAbstraction and takes only h, not (x0, h).
function make_grid_free(x0, h)
    uga = AB.UniformGridAbstraction

    if isdefined(uga, :MP) && isdefined(getfield(uga, :MP), :GridFree)
        return getfield(getfield(uga, :MP), :GridFree)(h)
    end

    error("Could not construct GridFree: MP.GridFree was not found in AB.UniformGridAbstraction.")
end


include("LTL_utils.jl")



############################################
################ Load LLM ##################
############################################
# Claude is queried lazily at runtime from `nl_to_spot_formula`.
# Set your Anthropic key before running this file, e.g. in the shell:
#   export ANTHROPIC_API_KEY="sk-ant-..."
# Optionally override the model with:
#   export CLAUDE_MODEL="claude-sonnet-4-5"
using HTTP
using JSON3

function load_env_file!(path::AbstractString = joinpath(@__DIR__, ".env"))
    isfile(path) || return false

    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "#") && continue
        occursin("=", line) || continue

        key, value = split(line, "="; limit = 2)
        key = strip(key)
        value = strip(value)

        # Remove optional surrounding quotes.
        if length(value) >= 2 && ((startswith(value, "\"") && endswith(value, "\"")) || (startswith(value, "'") && endswith(value, "'")))
            value = value[2:end-1]
        end

        # Do not overwrite values already exported in the shell.
        if !haskey(ENV, key)
            ENV[key] = value
        end
    end
    return true
end

load_env_file!()

const CLAUDE_API_URL = "https://api.anthropic.com/v1/messages"
const CLAUDE_API_KEY = get(ENV, "ANTHROPIC_API_KEY", "")
const CLAUDE_MODEL = get(ENV, "CLAUDE_MODEL", "claude-opus-4-7")
const CLAUDE_MAX_TOKENS = parse(Int, get(ENV, "CLAUDE_MAX_TOKENS", "256"))


const CLAUDE_SYSTEM_PROMPT = """
You translate natural-language robot task specifications into Linear Temporal Logic (LTL).
Return only one LTL formula and no explanation. Do not wrap the formula in Markdown code fences or backticks.
Use only these atomic propositions when applicable: blue, green, purple, brown, yellow, obs.
Use Spot-compatible syntax:
- G(prop_1) for always/globally
- F(prop_1) for eventually/finally
- X(prop_1) for next
- prop_1 U prop_2 for until
- !prop_1 for negation
- & for conjunction
- | for disjunction
- -> for implication
- <-> for equivalence
Always preserve the user's AP names exactly. Do not rename APs to prop_1, prop_2, etc.
"""

const CLAUDE_BACKTRANSLATION_PROMPT = """
You explain Linear Temporal Logic (LTL) formulas to non-expert robot users.
Return only a concise natural-language backtranslation of the formula. Do not mention syntax details unless necessary.
Use the same atomic proposition names appearing in the formula, such as blue, green, purple, brown, yellow, obs.
Be faithful to the formula and do not add requirements that are not present.
"""

const CLAUDE_VALIDATION_PROMPT = """
You are an assistant helping a user validate a robot task before planning.
Given the user's original request, the generated LTL formula, and a natural-language backtranslation of that formula, write a short validation message for the end user.
The message must:
1. State what the system understood the requirement to mean.
2. Ask the user to confirm whether this is correct.
3. Tell the user that if it is not correct, they should rephrase the requirement.
Do not use technical jargon. Do not mention Büchi automata, quotient systems, symbolic control, or controller synthesis.
Keep the message short.
"""


function call_claude_text(system_prompt::AbstractString, user_content::AbstractString; max_tokens::Int = CLAUDE_MAX_TOKENS)
    isempty(CLAUDE_API_KEY) && error("Missing Anthropic API key. Set ENV[\"ANTHROPIC_API_KEY\"] before running this file.")

    payload = Dict(
        "model" => CLAUDE_MODEL,
        "max_tokens" => max_tokens,
        #"temperature" => 0,
        "system" => String(system_prompt),
        "messages" => [
            Dict(
                "role" => "user",
                "content" => String(user_content),
            ),
        ],
    )

    headers = [
        "x-api-key" => CLAUDE_API_KEY,
        "anthropic-version" => "2023-06-01",
        "content-type" => "application/json",
    ]

    response = HTTP.post(CLAUDE_API_URL, headers, JSON3.write(payload); retry = false)
    if response.status < 200 || response.status >= 300
        error("Claude API request failed with status $(response.status): $(String(response.body))")
    end

    body = JSON3.read(String(response.body))

    if !haskey(body, :content) || isempty(body.content)
        error("Claude response did not contain any content: $(String(response.body))")
    end

    chunks = String[]
    for block in body.content
        if haskey(block, :type) && String(block.type) == "text" && haskey(block, :text)
            push!(chunks, String(block.text))
        end
    end

    isempty(chunks) && error("Claude response contained no text block: $(String(response.body))")
    return strip(join(chunks, "\n"))
end

function call_claude_ltl(nl_sentence::AbstractString)
    return call_claude_text(
        CLAUDE_SYSTEM_PROMPT,
        "Translate this natural-language task into one Spot-compatible LTL formula:\n$(String(nl_sentence))";
        max_tokens = CLAUDE_MAX_TOKENS,
    )
end

function call_claude_ltl_backtranslation(formula::AbstractString)
    return call_claude_text(
        CLAUDE_BACKTRANSLATION_PROMPT,
        "Backtranslate this LTL formula into concise natural language for a robot user:\n$(String(formula))";
        max_tokens = 256,
    )
end

function call_claude_validation_message(original_nl::AbstractString, formula::AbstractString, backtranslation::AbstractString)
    return call_claude_text(
        CLAUDE_VALIDATION_PROMPT,
        "Original user request:\n$(String(original_nl))\n\nGenerated LTL formula:\n$(String(formula))\n\nBacktranslation of the formula:\n$(String(backtranslation))";
        max_tokens = 256,
    )
end

function fallback_ltl_backtranslation(formula::AbstractString)
    ϕ = strip(String(formula))
    isempty(ϕ) && return "The system did not produce a task formula."
    return "The system interpreted the task as: $(ϕ)."
end

function fallback_validation_message(original_nl::AbstractString, formula::AbstractString, backtranslation::AbstractString)
    return "I understood your request as: $(String(backtranslation)) Please confirm whether this is correct. If not, please rephrase the requirement."
end

function build_requirement_validation_message(original_nl::AbstractString, formula::AbstractString)
    backtranslation = ""
    validation_message = ""

    try
        backtranslation = call_claude_ltl_backtranslation(formula)
    catch err
        @warn "LTL backtranslation LLM call failed; using fallback text." error=sprint(showerror, err)
        backtranslation = fallback_ltl_backtranslation(formula)
    end

    try
        validation_message = call_claude_validation_message(original_nl, formula, backtranslation)
    catch err
        @warn "Requirement validation LLM call failed; using fallback text." error=sprint(showerror, err)
        validation_message = fallback_validation_message(original_nl, formula, backtranslation)
    end

    return backtranslation, validation_message
end

function print_requirement_validation_summary(original_nl::AbstractString, formula::AbstractString, backtranslation::AbstractString, validation_message::AbstractString)
    println("\nRequirement validation draft")
    println("----------------------------")
    println("Original request: ", String(original_nl))
    println("Generated LTL:     ", String(formula))
    println("Backtranslation:   ", String(backtranslation))
    println("User message:      ", String(validation_message))
    println("----------------------------\n")
    return nothing
end

# Global constant for the input box (used by input grid helpers)
const U_BOX = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))


# Global time step used consistently in abstraction, cost, and simulation

const DT = 0.30

# Safety margin used to enlarge obstacles during abstraction and labeling.
# This keeps the synthesized/executed trajectory away from obstacle boundaries.
const OBSTACLE_MARGIN = 0.12
const DEFAULT_QUOTIENT_SAVE_PATH = joinpath(@__DIR__, "saved_ap_label_quotient.jls")

if !isdefined(Main, :LTL2Automata)
    include("src/LTL/LTL2Automata.jl")
end

import .LTL2Automata
const L2A = LTL2Automata

using .LTL2Automata: TSFromAdj, letter_of, nextstate, AbstractTS, translate_ltl_buchi,
                    Automaton, is_accepting, ts_from_symbolic_onepass,
                    synthesize_buchi_controller, simulate_buchi_policy,
                    precompute_labels_and_forbidden,
                    APLabeler, labeler_and_inner_sets
# Toggle to disable heavy automaton diagnostics/printing
const AUT_DIAGNOSTICS = false





########### Interactive mode: no startup LTL specification ###########
######################################################################
# This file should build the abstraction / TS / labels first, then wait
# for the user to provide a natural-language specification at runtime.
# NL -> LTL translation will happen later in the interactive layer.
const DEFAULT_STARTUP_LTL = ""
######################################################################




##########################################################
########### Define system and buil abstraction ###########
##########################################################

include("System_helper.jl")



x_start = [4.8, 5.0 , 0.0]
x0 = SVector(0.0, 0.0, pi/2);
hx = SVector(0.20, 0.20, 0.3);
u0 = SVector(0.0, 0.0);
hu = SVector(0.3, 0.3);
start = SVector(0.0); # SVector(1.0);

include("abstraction_helper.jl")

optimizer = MOI.instantiate(AB.UniformGridAbstraction.Optimizer)
concrete_problem = problem(x1_lb, x1_ub, x2_lb, x2_ub)
MOI.set(optimizer, MOI.RawOptimizerAttribute("concrete_problem"), concrete_problem)
MOI.set(optimizer, MOI.RawOptimizerAttribute("state_grid"), make_grid_free(x0, hx))
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("input_grid"),
    make_grid_free(u0, hu),
)
MOI.set(
    optimizer,  
    MOI.RawOptimizerAttribute("jacobian_bound"),
    jacobian_bound(),
)
MOI.set(optimizer, MOI.RawOptimizerAttribute("time_step"), DT)
MOI.set(
    optimizer,
    MOI.RawOptimizerAttribute("approx_mode"),
    AB.UniformGridAbstraction.GROWTH,
)

println("Starting abstraction…")
t = @elapsed begin
    MOI.optimize!(optimizer)
end
println("Abstraction time: $t seconds")

abstract_system = MOI.get(optimizer, MOI.RawOptimizerAttribute("abstract_system"));

########################################
########### Label Each State ###########
########################################
# This section is to Assign an Atomic Proposition to each region

# Define the regions to label
X = concrete_problem.region
θ_min = X.lb[3]
θ_max = X.ub[3]
blue = UT.HyperRectangle(
    SVector(6.0, 7.5, θ_min),
    SVector(9.0, 9.5, θ_max),
)

green = UT.HyperRectangle(
    SVector(5.0, 1.0, θ_min),
    SVector(7.5, 2.5, θ_max),
)

purple = UT.HyperRectangle(
    SVector(7.0, 3.5, θ_min),
    SVector(8.5, 6.5, θ_max),
)
brown = UT.HyperRectangle(
    SVector(4.0, 7.0, θ_min),
    SVector(7.5, 8.0, θ_max),
)
yellow = UT.HyperRectangle(
    SVector(0.5, 1.0, θ_min),
    SVector(2.5, 3.5, θ_max),
)





obstacle_rects = UT.HyperRectangle[]
for i in eachindex(x1_lb)
    push!(obstacle_rects,
        UT.HyperRectangle(
            SVector(
                max(X.lb[1], x1_lb[i] - OBSTACLE_MARGIN),
                max(X.lb[2], x2_lb[i] - OBSTACLE_MARGIN),
                θ_min,
            ),
            SVector(
                min(X.ub[1], x1_ub[i] + OBSTACLE_MARGIN),
                min(X.ub[2], x2_ub[i] + OBSTACLE_MARGIN),
                θ_max,
            ),
        )
    )
end

regions = Dict(
    :obs  => obstacle_rects,
    :blue => blue,
    :green  => green,
    :purple => purple,
    :brown  => brown,
    :yellow => yellow,
)


# which APs are forbidden (used in G(!obs), etc.)
forbidden_aps = [:obs]

# build OUTER-based labeler + INNER sets for all non-forbidden APs
Label, inner_sets = labeler_and_inner_sets(
    abstract_system;
    regions       = regions,
    forbidden_aps = forbidden_aps,
)


###############################################
########### Build Transition System ###########
###############################################


#=
function _post_automaton(abstract_system::SY.SymbolicModelList,x::Int, u::Int)
    autom = abstract_system.autom  # this is an AutomatonList
    targets = Int[]
    # transitions are (target, source, symbol); the "tail" we fix on is (source, symbol)
    UT.fix_and_eliminate_tail!(targets, autom.transitions, (x, u))
    return targets
end
=#

x0_cont = SVector(4.8, 5.0 , pi/2)
eps = SVector(0.10, 0.10, 0.10)
X_start = UT.HyperRectangle(x0_cont - eps, x0_cont + eps)

# Build TS adjacency once from the symbolic model
# (ts_from_symbolic_onepass in LTL2Automata does not accept a u_of keyword)
_ts_tmp = ts_from_symbolic_onepass(
    abstract_system,
    X_start;
    reachable_only = false,
    x0_abs = nothing,
)

# Globally remove obstacle states from the transition system.
# This ensures that the controller never plans through black obstacle regions,
# even when the user requirement does not explicitly mention `obs`.
function prune_forbidden_successors(ts_local::TSFromAdj, Label, forbidden_aps::Vector{Symbol})
    nx_local = L2A.nstates(ts_local)
    nu_local = L2A.ninputs(ts_local)

    is_forbidden_state(x::Int) = any(ap -> ap in Label(x), forbidden_aps)

    safe_adj = [
        [
            filter(xp -> !is_forbidden_state(xp), L2A.post(ts_local, x, u))
            for u in 1:nu_local
        ]
        for x in 1:nx_local
    ]

    n_removed = 0
    for x in 1:nx_local
        for u in 1:nu_local
            n_removed += length(L2A.post(ts_local, x, u)) - length(safe_adj[x][u])
        end
    end

    println("Pruned ", n_removed, " transitions entering forbidden obstacle states.")
    return L2A.TSFromAdj(safe_adj, L2A.initset(ts_local))
end

_ts_tmp = prune_forbidden_successors(_ts_tmp, Label, forbidden_aps)

# Backward-compat alias (older scripts used ts_from_symbolic)
const ts_from_symbolic = ts_from_symbolic_onepass

# Force X0 to be the *single* abstract state whose cell contains x0_cont.
# We query a tiny box of half-cell size around x0_cont and use DO.CENTER so
# the returned state corresponds to the representative point inside that cell.
X0_cell = UT.HyperRectangle(x0_cont .- (hx ./ 2), x0_cont .+ (hx ./ 2)) 
X0_center_states = SY.get_states_from_set(abstract_system, X0_cell, AB.UniformGridAbstraction.MP.OUTER)

if isempty(X0_center_states)
    # Fallback: OUTER query (should almost always return something)
    X0_center_states = SY.get_states_from_set(abstract_system, X0_cell, AB.UniformGridAbstraction.MP.OUTER)
end

if isempty(X0_center_states)
    error("Could not map x0_cont to any abstract state. Check that x0_cont is inside X and hx is consistent.")
end

x0_abs_start = X0_center_states[1]
println("Forced TS initial state from x0_cont: x0_abs_start=", x0_abs_start, ", label=", Label(x0_abs_start))


# Final TS uses the same adjacency but a single initial state at x0_cont

ts = L2A.TSFromAdj(_ts_tmp.Adj, [x0_abs_start])
if any(ap -> ap in Label(x0_abs_start), forbidden_aps)
    error("The initial state is labeled as forbidden/obstacle. Move x0_cont or reduce OBSTACLE_MARGIN.")
end

# Precompute TS labels once at startup so the AP-label quotient can be built and saved immediately.
const TASK_APS = Symbol[:obs, :blue, :green, :purple, :brown, :yellow]
const LABELS_TS_STARTUP, _ = precompute_labels_and_forbidden(ts, Label, TASK_APS)


# =============================================================================
# Phase 1: layered architecture
#   - backend: expensive objects built once
#   - translation: NL -> Spot-LTL
#   - Büchi build: Spot-LTL -> automaton + TS labels
# =============================================================================


Base.@kwdef mutable struct PlannerBackend
    claude_model::String
    claude_api_key::String
    claude_api_url::String
    claude_max_tokens::Int
    abstract_system
    ts::TSFromAdj
    Label
    inner_sets::Dict{Symbol, Set{Int}}
    x0_cont::SVector{3,Float64}
    x0_abs_start::Int
    hx::SVector{3,Float64}
end


Base.@kwdef mutable struct BuchiBuildResult
    nl::String
    formula::String
    aut::Automaton
    labels_ts::Vector{Set{Symbol}}
    backtranslation::String = ""
    validation_message::String = ""
end


Base.@kwdef mutable struct QuotientCheckResult
    compatible::Bool
    message::String
    accepting_product_state::Union{Nothing, Tuple{Int,Int}} = nothing
    nquotient_states::Int = 0
    nproduct_states_reached::Int = 0
end


"""Build a coarse quotient by grouping concrete TS states with exactly the same AP label set.

This is a lightweight compatibility abstraction. It is intentionally existential:
a quotient edge q -> q' exists if at least one concrete transition between the
corresponding AP-label classes exists under some input.
"""
function build_ap_label_quotient(ts_local::TSFromAdj, labels_ts::Vector{Set{Symbol}})
    nx = L2A.nstates(ts_local)
    nu = L2A.ninputs(ts_local)

    label_keys = [Tuple(sort!(collect(labels_ts[x]))) for x in 1:nx]
    key_to_q = Dict{Tuple{Vararg{Symbol}}, Int}()
    state_to_q = Vector{Int}(undef, nx)
    quotient_labels = Vector{Set{Symbol}}()

    for x in 1:nx
        key = label_keys[x]
        q = get(key_to_q, key, 0)
        if q == 0
            q = length(quotient_labels) + 1
            key_to_q[key] = q
            push!(quotient_labels, Set{Symbol}(key))
        end
        state_to_q[x] = q
    end

    nq = length(quotient_labels)
    adj_sets = [Set{Int}() for _ in 1:nq]

    for x in 1:nx
        qx = state_to_q[x]
        for u in 1:nu
            for xp in L2A.post(ts_local, x, u)
                push!(adj_sets[qx], state_to_q[xp])
            end
        end
    end

    quotient_adj = [sort!(collect(succs)) for succs in adj_sets]
    quotient_x0 = unique!(sort!([state_to_q[x0] for x0 in L2A.initset(ts_local)]))

    return quotient_adj, quotient_x0, quotient_labels, state_to_q, key_to_q
end

function save_ap_label_quotient(
    ts_local::TSFromAdj,
    labels_ts::Vector{Set{Symbol}};
    output_path::String = DEFAULT_QUOTIENT_SAVE_PATH,
)
    quotient_adj, quotient_x0, quotient_labels, state_to_q, key_to_q = build_ap_label_quotient(ts_local, labels_ts)

    payload = Dict(
        "quotient_adj" => quotient_adj,
        "quotient_x0" => quotient_x0,
        "quotient_labels" => quotient_labels,
        "state_to_q" => state_to_q,
        "label_key_to_q" => Dict(string(k) => v for (k, v) in key_to_q),
        "num_quotient_states" => length(quotient_labels),
        "num_ts_states" => L2A.nstates(ts_local),
        "num_ts_inputs" => L2A.ninputs(ts_local),
    )

    open(output_path, "w") do io
        serialize(io, payload)
    end

    println("Saved AP-label quotient to: ", output_path)
    return output_path
end


function load_ap_label_quotient(path::String = DEFAULT_QUOTIENT_SAVE_PATH)
    isfile(path) || error("Saved quotient file not found: $(path)")
    open(path, "r") do io
        return deserialize(io)
    end
end

save_ap_label_quotient(ts, LABELS_TS_STARTUP)


function _automaton_ap_names(aut_local::Automaton)
    for fname in (:ap, :aps, :AP, :atomic_propositions)
        if hasfield(typeof(aut_local), fname)
            return Vector{String}(getfield(aut_local, fname))
        end
    end
    error("Could not find the atomic-proposition list field in Automaton. Expected one of: :ap, :aps, :AP, :atomic_propositions.")
end

function _buchi_successor_after_quotient_label(aut_local::Automaton, qaut::Int, label::Set{Symbol})
    # LTL2Automata.letter_of has signature letter_of(::Set{Symbol}, ::Vector{String}).
    # It converts the AP label set of a quotient state into the letter expected
    # by the Büchi transition function.
    letter = letter_of(label, _automaton_ap_names(aut_local))
    return nextstate(aut_local, qaut, letter)
end

"""Check whether the synchronous product of the AP-label quotient and Büchi automaton has an accepting run.

The check computes the reachable product graph and then searches for a reachable
SCC containing an accepting Büchi state and a cycle. This is a fast compatibility
precheck, not a proof of robust controllability on the full abstraction.
"""
function check_formula_on_ap_quotient(ts_local::TSFromAdj, aut_local::Automaton, labels_ts::Vector{Set{Symbol}})
    quotient_adj, quotient_x0, quotient_labels, _, _ = build_ap_label_quotient(ts_local, labels_ts)
    nq = length(quotient_labels)
    isempty(quotient_x0) && return QuotientCheckResult(
        compatible = false,
        message = "Quotient compatibility check failed: empty quotient initial set.",
        nquotient_states = nq,
        nproduct_states_reached = 0,
    )

    # Use reached product pairs as graph nodes. Product transition:
    #   (q, b) -> (q_next, b_next), where b_next reads label(q_next).
    start_pairs = Tuple{Int,Int}[]
    for q0 in quotient_x0
        b0_after_label = _buchi_successor_after_quotient_label(aut_local, aut_local.init, quotient_labels[q0])
        b0_after_label === nothing && continue
        push!(start_pairs, (q0, b0_after_label))
    end

    isempty(start_pairs) && return QuotientCheckResult(
        compatible = false,
        message = "Quotient compatibility check failed: no valid Büchi transition from the initial quotient label.",
        nquotient_states = nq,
        nproduct_states_reached = 0,
    )

    pair_to_id = Dict{Tuple{Int,Int}, Int}()
    id_to_pair = Tuple{Int,Int}[]
    product_adj = Vector{Vector{Int}}()
    queue = Tuple{Int,Int}[]

    function ensure_pair!(pair::Tuple{Int,Int})
        id = get(pair_to_id, pair, 0)
        if id == 0
            id = length(id_to_pair) + 1
            pair_to_id[pair] = id
            push!(id_to_pair, pair)
            push!(product_adj, Int[])
            push!(queue, pair)
        end
        return id
    end

    for pair in start_pairs
        ensure_pair!(pair)
    end

    head = 1
    while head <= length(queue)
        pair = queue[head]
        head += 1
        src_id = pair_to_id[pair]
        q, b = pair
        for qnext in quotient_adj[q]
            bnext = _buchi_successor_after_quotient_label(aut_local, b, quotient_labels[qnext])
            bnext === nothing && continue
            dst = (qnext, bnext)
            dst_id = ensure_pair!(dst)
            push!(product_adj[src_id], dst_id)
        end
        unique!(product_adj[src_id])
    end

    nprod = length(id_to_pair)
    nprod == 0 && return QuotientCheckResult(
        compatible = false,
        message = "Quotient compatibility check failed: no reachable product states.",
        nquotient_states = nq,
        nproduct_states_reached = 0,
    )

    # Tarjan SCC over reachable product graph.
    index = 0
    indices = fill(0, nprod)
    lowlink = fill(0, nprod)
    onstack = falses(nprod)
    stack = Int[]
    compatible_pair = Ref{Union{Nothing, Tuple{Int,Int}}}(nothing)

    function strongconnect(v::Int)
        index += 1
        indices[v] = index
        lowlink[v] = index
        push!(stack, v)
        onstack[v] = true

        for w in product_adj[v]
            if indices[w] == 0
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elseif onstack[w]
                lowlink[v] = min(lowlink[v], indices[w])
            end
        end

        if lowlink[v] == indices[v]
            scc = Int[]
            while true
                w = pop!(stack)
                onstack[w] = false
                push!(scc, w)
                w == v && break
            end

            has_accepting = any(id -> is_accepting(aut_local, id_to_pair[id][2]), scc)
            has_cycle = length(scc) > 1 || any(w -> w == scc[1], product_adj[scc[1]])
            if has_accepting && has_cycle && compatible_pair[] === nothing
                accepting_id = first(filter(id -> is_accepting(aut_local, id_to_pair[id][2]), scc))
                compatible_pair[] = id_to_pair[accepting_id]
            end
        end
        return nothing
    end

    for v in 1:nprod
        compatible_pair[] !== nothing && break
        indices[v] == 0 && strongconnect(v)
    end

    if compatible_pair[] === nothing
        return QuotientCheckResult(
            compatible = false,
            message = "Requirement is not compatible with the coarse AP-label quotient: no reachable accepting Büchi cycle exists in the quotient product.",
            accepting_product_state = nothing,
            nquotient_states = nq,
            nproduct_states_reached = nprod,
        )
    end

    return QuotientCheckResult(
        compatible = true,
        message = "Requirement passed the coarse AP-label quotient compatibility check.",
        accepting_product_state = compatible_pair[],
        nquotient_states = nq,
        nproduct_states_reached = nprod,
    )
end

function assert_quotient_compatible!(backend::PlannerBackend, build::BuchiBuildResult)
    result = check_formula_on_ap_quotient(backend.ts, build.aut, build.labels_ts)
    println("Quotient compatibility check: ", result.message)
    println("  quotient states: ", result.nquotient_states)
    println("  reachable product states: ", result.nproduct_states_reached)
    if !result.compatible
        error(result.message)
    end
    return result
end



"""Collect all expensive objects built at startup into one backend object."""
function build_backend_from_globals()
    return PlannerBackend(
        claude_model = CLAUDE_MODEL,
        claude_api_key = CLAUDE_API_KEY,
        claude_api_url = CLAUDE_API_URL,
        claude_max_tokens = CLAUDE_MAX_TOKENS,
        abstract_system = abstract_system,
        ts = ts,
        Label = Label,
        inner_sets = inner_sets,
        x0_cont = SVector{3,Float64}(x0_cont...),
        x0_abs_start = x0_abs_start,
        hx = SVector{3,Float64}(hx...),
    )
end

"""Translate a user-provided natural-language requirement into Spot-style LTL.
This keeps the LLM loaded at startup, but defers actual translation until the
user enters a requirement at runtime.
"""
function nl_to_spot_formula(backend::PlannerBackend, nl_sentence::AbstractString)
    nl_sentence = strip(String(nl_sentence))
    isempty(nl_sentence) && error("Empty NL specification.")

    println("NL (raw):   ", nl_sentence)

    # Keep the user's atomic proposition names unchanged.
    # Claude is expected to produce LTL directly over the environment AP names,
    # e.g., blue, green, purple, brown, yellow, obs.
    decoded = call_claude_ltl(nl_sentence)
    println("Claude output (raw LTL-like): ", decoded)

    ϕ_local = sanitize_spot_formula_candidate(String(strip(llm_to_spot_ltl(decoded))))

    # Minimal safety net only for clearly broken outputs such as
    # empty strings, punctuation-only junk, or incomplete formulas like `X`.
    # Valid formulas are left untouched and parsed by Spot directly.
    if is_obviously_incomplete_spot_formula(ϕ_local) || is_obviously_invalid_spot_formula(ϕ_local)
        @warn "Claude produced a broken Spot formula; using fallback NL parser." nl_sentence decoded ϕ_local
        ϕ_local = fallback_nl_to_spot_formula(nl_sentence)
    end

    println("Spot-style ϕ from LLM: ", ϕ_local)
    return ϕ_local
end

# Backward-compatible convenience wrapper
function nl_to_spot_formula(nl_sentence::AbstractString)
    return nl_to_spot_formula(BACKEND, nl_sentence)
end

"""Build the Büchi-side objects for one NL requirement using an already-built backend."""
function build_buchi_from_nl(backend::PlannerBackend, nl_sentence::AbstractString)
    ϕ_local = String(nl_to_spot_formula(backend, nl_sentence))
    isempty(strip(ϕ_local)) && error("Refusing to build Büchi automaton from an empty formula.")

    backtranslation, validation_message = build_requirement_validation_message(nl_sentence, ϕ_local)
    print_requirement_validation_summary(nl_sentence, ϕ_local, backtranslation, validation_message)

    println("Translating LTL to Büchi automaton…")
    _t_aut = @elapsed begin
        aut_local = L2A.translate_ltl_buchi(ϕ_local)
        labels_local, _ = precompute_labels_and_forbidden(backend.ts, backend.Label, Symbol[])
        global LAST_BUCHI_BUILD = BuchiBuildResult(
            nl = String(nl_sentence),
            formula = ϕ_local,
            aut = aut_local,
            labels_ts = labels_local,
            backtranslation = String(backtranslation),
            validation_message = String(validation_message),
        )
        global LAST_QUOTIENT_CHECK = assert_quotient_compatible!(backend, LAST_BUCHI_BUILD)
    end
    println("LTL→Büchi + labels time: ", _t_aut, " seconds")
    return LAST_BUCHI_BUILD
end



using Plots

X = concrete_problem.region    # HyperRectangle
L = X.ub .- X.lb               # side lengths

nx1 = Int(round(L[1] / hx[1]))   # number of cells along x1
nx2 = Int(round(L[2] / hx[2]))
nx3 = Int(round(L[3] / hx[3]))
nx  = (nx1, nx2, nx3)




# -----------------------------------------------------------------------------
# Fast center lookup for plotting
# -----------------------------------------------------------------------------
# Newer Dionysos versions do not expose a per-state “cell” accessor.
# However, SymbolicModelList stores a bijection between abstract state ids
# and grid positions via `xint2pos` / `xpos2int`, and Xdom stores the grid.
# We use that to recover the (x,y) center of any abstract state.


#
# Extract (x0, hx) from a Dionysos grid object (GridFree).
# Tries common names; if that fails, scans fields for StaticVectors/Tuples.
#






function plot_abstract_cells!(fig, state_map)
    first = true
    for rect in values(state_map)
        xs = [rect.lb[1], rect.ub[1], rect.ub[1], rect.lb[1]]
        ys = [rect.lb[2], rect.lb[2], rect.ub[2], rect.ub[2]]
        plot!(fig, Shape(xs, ys);
              fillcolor = :lightgray,
              fillalpha = 0.1,
              linecolor = :lightgray,
              linealpha = 0.3,
              label = first ? "abstract cells" : "")
        first = false
    end
end

function plot_environment()
    fig = plot(; aspect_ratio = :equal, legend = :bottomleft)

    # Concrete domain
    plot!(fig, concrete_problem.region; color = :grey, opacity = 0.5, label = "")
    # Draw abstract cells directly from the abstraction domain. This is much faster
    # than plotting each cell individually and gives the familiar cyan grid.
    plot!(fig, abstract_system.Xdom;
          color = :deepskyblue,
          opacity = 0.75,
          efficient = false,
          linewidth = 0.7,
          label = false)

    # Regions you already defined
    blue_shape = Shape(
        [6.0, 9.0, 9.0, 6.0],
        [7.5, 7.5, 9.5, 9.5],
    )
    green_shape = Shape(
        [5.0, 7.5, 7.5, 5.0],
        [1.0, 1.0, 2.5, 2.5],
    )
    purple_shape = Shape(
        [7.0, 8.5, 8.5, 7.0],
        [3.5, 3.5, 6.5, 6.5],
    )
    brown_shape = Shape(
        [4.0, 7.5, 7.5, 4.0],
        [7.0, 7.0, 8.0, 8.0],
    )
    yellow_shape = Shape(
        [0.5, 2.5, 2.5, 0.5],
        [1.0, 1.0, 3.5, 3.5],
    )
    plot!(fig, blue_shape;  fillcolor=:blue,  fillalpha=0.7, linecolor=:black, label=false)
    plot!(fig, green_shape; fillcolor=:green, fillalpha=0.7, linecolor=:black, label=false)
    plot!(fig, purple_shape; fillcolor=:purple, fillalpha=0.7, linecolor=:black, label=false)
    plot!(fig, brown_shape; fillcolor=:brown, fillalpha=0.7, linecolor=:black, label=false)
    plot!(fig, yellow_shape; fillcolor=:yellow, fillalpha=0.7, linecolor=:black, label=false)
    # Obstacles (projection to x1,x2)
    for i in eachindex(x1_lb)
        obs_shape = Shape(
            [x1_lb[i], x1_ub[i], x1_ub[i], x1_lb[i]],
            [x2_lb[i], x2_lb[i], x2_ub[i], x2_ub[i]],
        )
        plot!(fig, obs_shape;
              fillcolor = :black,
              fillalpha = 0.9,
              linecolor = :black,
              linewidth = 1.8,
              label = (i == 1 ? "obstacle" : ""))
    end

    return fig
end



# ------------------------------------------------------------------------
# Helper: map input index to the actual abstract input vector, and discrete plan follower
# ------------------------------------------------------------------------

"""
    input_from_index(u_idx, abstract_system)

Return the *actual* 2D control vector corresponding to the abstract
input index `u_idx`, using `abstract_system.Udom`, which is a
`CustomList{2,Float64}` of `SVector{2,Float64}`. This avoids guessing
a grid layout and uses exactly the same inputs as the abstraction.
"""
 

# --------------------------------------------------------------------------
# Wait for a natural-language requirement from the user before any LTL/Büchi work.
# In VS Code / REPL setups, reading explicitly from `stdin` is more reliable
# than a bare `readline()` at top level.
# --------------------------------------------------------------------------

function wait_for_user_nl_spec()
    while true
        println("\nEnter a natural-language requirement for the planner:")
        print("NL spec> ")
        flush(stdout)

        line = strip(readline(stdin))
        !isempty(line) && return line

        println("Empty input received. Please type a natural-language requirement.")
    end
end

function run_buchi_from_user_nl!()
    global USER_NL_SPEC = wait_for_user_nl_spec()
    result = build_buchi_from_nl(BACKEND, USER_NL_SPEC)

    global ϕ = result.formula
    global aut = result.aut
    global labels_ts = result.labels_ts
    global LAST_BUCHI_BUILD = result

    if !(@isdefined(RUNTIME)) || RUNTIME === nothing
        global RUNTIME = init_runtime(BACKEND)
    end
    reset_runtime_plan!(RUNTIME)
    install_buchi_build!(RUNTIME, result)

    println("Built Büchi objects for NL requirement.")
    println("Formula: ", ϕ)
    println("Runtime updated. Inspect with: describe_runtime(RUNTIME)")
    return nothing
end

global BACKEND = build_backend_from_globals()
println("Backend object created.")
println("You can inspect it with: BACKEND")

# =============================================================================
# Phase 2: persistent runtime object
#   Stores the current interactive session state independently from the backend.
# =============================================================================


if !isdefined(Main, :PlannerRuntime)
    Base.@kwdef mutable struct PlannerRuntime
        backend::PlannerBackend
        current_abs::Int
        current_state::SVector{3,Float64}
        current_formula::String = ""
        current_aut = nothing
        current_labels_ts::Union{Nothing, Vector{Set{Symbol}}} = nothing
        current_plan = nothing
        plan_token::Int = 0
        latest_request_nl::String = ""
        active_task = nothing
        is_executing::Bool = false
        trajectory_abs::Vector{Int} = Int[]
        trajectory_xy::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]
    end
end

# Backward-compatible side storage for runtime metadata when an older
# PlannerRuntime type is already loaded in the REPL and cannot be redefined.
if !isdefined(Main, :_RUNTIME_LATEST_REQUEST)
    const _RUNTIME_LATEST_REQUEST = IdDict{Any,String}()
end
if !isdefined(Main, :_RUNTIME_ACTIVE_TASK)
    const _RUNTIME_ACTIVE_TASK = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_PLAN_STATUS)
    const _RUNTIME_PLAN_STATUS = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_PLAN_XS)
    const _RUNTIME_PLAN_XS = IdDict{Any,Vector{Int}}()
end
if !isdefined(Main, :_RUNTIME_PLAN_IQS)
    const _RUNTIME_PLAN_IQS = IdDict{Any,Vector{Int}}()
end
if !isdefined(Main, :_RUNTIME_PLAN_US)
    const _RUNTIME_PLAN_US = IdDict{Any,Vector{Int}}()
end
if !isdefined(Main, :_RUNTIME_WINP)
    const _RUNTIME_WINP = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_POLP)
    const _RUNTIME_POLP = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_NEXTQP)
    const _RUNTIME_NEXTQP = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_INFO)
    const _RUNTIME_INFO = IdDict{Any,Any}()
end
if !isdefined(Main, :_RUNTIME_EXEC_TASK)
    const _RUNTIME_EXEC_TASK = IdDict{Any,Any}()
end
if !isdefined(Main, :_UI_PATH_OBS)
    const _UI_PATH_OBS = IdDict{Any,Any}()
end

"""Initialize the interactive runtime from an already-built backend."""
function init_runtime(backend::PlannerBackend)
    x_abs = backend.x0_abs_start
    x_cont = backend.x0_cont
    x1, x2 = center2d(x_abs)

    rt = PlannerRuntime(
        backend = backend,
        current_abs = x_abs,
        current_state = x_cont,
        trajectory_abs = [x_abs],
        trajectory_xy = [(x1, x2)],
    )
    return rt
end


"""Reset runtime bookkeeping when a new plan is about to be installed."""
function reset_runtime_plan!(rt::PlannerRuntime)
    rt.current_formula = ""
    rt.current_aut = nothing
    rt.current_labels_ts = nothing
    rt.current_plan = nothing
    rt.is_executing = false
    clear_runtime_planning!(rt)

    # Start a fresh displayed trajectory from the current state.
    x1, x2 = center2d(rt.current_abs)
    rt.trajectory_abs = [rt.current_abs]
    rt.trajectory_xy = [(Float64(x1), Float64(x2))]
    return rt
end

"""Start a new replanning request and invalidate all older pending results."""
function begin_replan!(rt::PlannerRuntime, nl_text::AbstractString)
    rt.plan_token += 1
    _RUNTIME_LATEST_REQUEST[rt] = String(nl_text)
    rt.is_executing = false
    rt.current_plan = nothing
    return rt.plan_token
end

"""Return true iff `token` is still the latest request token."""
is_current_plan_token(rt::PlannerRuntime, token::Int) = (token == rt.plan_token)

"""Install the result of one NL->LTL->Büchi build into the runtime."""
function install_buchi_build!(rt::PlannerRuntime, build::BuchiBuildResult)
    rt.current_formula = build.formula
    rt.current_aut = build.aut
    rt.current_labels_ts = build.labels_ts
    rt.current_plan = build
    return rt
end

"""Clear any previously stored synthesized plan for this runtime."""
function clear_runtime_planning!(rt::PlannerRuntime)
    delete!(_RUNTIME_PLAN_STATUS, rt)
    delete!(_RUNTIME_PLAN_XS, rt)
    delete!(_RUNTIME_PLAN_IQS, rt)
    delete!(_RUNTIME_PLAN_US, rt)
    delete!(_RUNTIME_WINP, rt)
    delete!(_RUNTIME_POLP, rt)
    delete!(_RUNTIME_NEXTQP, rt)
    delete!(_RUNTIME_INFO, rt)
    return rt
end

"""Synthesize and simulate a Büchi plan from the runtime's current abstract state."""
function plan_from_runtime!(rt::PlannerRuntime; max_steps::Int = 1200, stop_at_accept::Bool = true, semantics = :robust)
    rt.current_aut === nothing && error("Runtime has no current automaton. Submit an NL requirement first.")
    rt.current_labels_ts === nothing && error("Runtime has no current TS labels. Submit an NL requirement first.")

    global LAST_QUOTIENT_CHECK = check_formula_on_ap_quotient(rt.backend.ts, rt.current_aut, rt.current_labels_ts)
    println("Quotient compatibility check: ", LAST_QUOTIENT_CHECK.message)
    println("  quotient states: ", LAST_QUOTIENT_CHECK.nquotient_states)
    println("  reachable product states: ", LAST_QUOTIENT_CHECK.nproduct_states_reached)
    if !LAST_QUOTIENT_CHECK.compatible
        error(LAST_QUOTIENT_CHECK.message)
    end

    ts_local = rt.backend.ts
    aut_local = rt.current_aut
    labels_local = rt.current_labels_ts
    x0_local = rt.current_abs
    q0_local = aut_local.init

    println("Planning from current abstract state ...")
    _t_plan = @elapsed begin
        WinP, polP, next_qP, infoP = synthesize_buchi_controller(ts_local, aut_local, labels_local; semantics = semantics)
        xs, iqs, us, status = simulate_buchi_policy(
            ts_local, aut_local, labels_local, polP, next_qP;
            x0 = x0_local,
            q0 = q0_local,
            max_steps = max_steps,
            stop_at_accept = stop_at_accept,
        )

        _RUNTIME_WINP[rt] = WinP
        _RUNTIME_POLP[rt] = polP
        _RUNTIME_NEXTQP[rt] = next_qP
        _RUNTIME_INFO[rt] = infoP
        _RUNTIME_PLAN_XS[rt] = xs
        _RUNTIME_PLAN_IQS[rt] = iqs
        _RUNTIME_PLAN_US[rt] = us
        _RUNTIME_PLAN_STATUS[rt] = status
    end
    println("Planning time: ", _t_plan, " seconds")
    return get(_RUNTIME_PLAN_STATUS, rt, nothing)
end

"""Print a short summary of the currently stored synthesized plan."""
function describe_plan(rt::PlannerRuntime)
    status = get(_RUNTIME_PLAN_STATUS, rt, nothing)
    xs = get(_RUNTIME_PLAN_XS, rt, Int[])
    us = get(_RUNTIME_PLAN_US, rt, Int[])
    println("PlannerRuntime plan summary:")
    println("  status        = ", status === nothing ? "<none>" : status)
    println("  nstates       = ", length(xs))
    println("  ninputs       = ", length(us))
    if !isempty(xs)
        println("  start_state   = ", first(xs))
        println("  end_state     = ", last(xs))
    end
    return nothing
end

"""Print a short summary of the current runtime state."""
function describe_runtime(rt::PlannerRuntime)
    println("PlannerRuntime summary:")
    println("  current_abs      = ", rt.current_abs)
    println("  current_state    = ", rt.current_state)
    println("  current_formula  = ", isempty(rt.current_formula) ? "<none>" : rt.current_formula)
    latest_req = get(_RUNTIME_LATEST_REQUEST, rt, "")
    println("  latest_request   = ", isempty(latest_req) ? "<none>" : latest_req)
    println("  has automaton?   = ", rt.current_aut !== nothing)
    println("  has labels?      = ", rt.current_labels_ts !== nothing)
    println("  has plan?        = ", rt.current_plan !== nothing)
    println("  plan_token       = ", rt.plan_token)
    println("  is_executing     = ", rt.is_executing)
    println("  has_exec_task    = ", haskey(_RUNTIME_EXEC_TASK, rt))
    println("  trajectory_len   = ", length(rt.trajectory_abs))
    plan_status = get(_RUNTIME_PLAN_STATUS, rt, nothing)
    println("  plan_status     = ", plan_status === nothing ? "<none>" : plan_status)
    return nothing
end

global RUNTIME = init_runtime(BACKEND)
println("Runtime object created.")
println("You can inspect it with: describe_runtime(RUNTIME)")
println("The AP-label quotient is saved automatically at startup. You can also use save_ap_label_quotient(ts, LABELS_TS_STARTUP) and load_ap_label_quotient().")

#println("Translating LTL to Büchi automaton…")
#_t_aut = @elapsed begin
#    global aut = L2A.translate_ltl_buchi(ϕ)
#end
#println("LTL→Büchi time: ", _t_aut, " seconds")


#println("Precomputing TS labels…")
#_t_labels = @elapsed begin
#    global labels_ts, _ = precompute_labels_and_forbidden(ts, Label, Symbol[])
#end
#println("Label precompute time: ", _t_labels, " seconds")

#
# --------------------------------------------------------------------------
# Static demo path (disabled by default).
# In interactive mode, nothing below should run until the user explicitly
# provides NL and triggers the Büchi pipeline.
# --------------------------------------------------------------------------
const RUN_BUCHI_GAME = false
const GAME_SEMANTICS = :robust
const RUN_STATIC_DEMO = false

if !RUN_STATIC_DEMO
    println("Static Büchi/TS demo is disabled at startup. Use run_buchi_from_user_nl!() after include(...).")
elseif RUN_BUCHI_GAME && @isdefined(labels_ts) && @isdefined(aut)
let
    println("Projected label nonempty states: ", count(!isempty, labels_ts), " / ", length(labels_ts))

    # Solve Büchi game on product
    _t_synth = @elapsed begin
        global WinP, polP, next_qP, infoP = synthesize_buchi_controller(ts, aut, labels_ts; semantics = :robust)
    end
    println("Büchi game synthesis time: ", _t_synth, " seconds")

    nx_ts = L2A.nstates(ts)
    nqP   = infoP.nq
    np    = nx_ts * nqP
    println("Büchi: winning product states (robust) = ", count(WinP), " / ", np)

    x0_abs = ts.X0[1]
    q0_aut = aut.init

    # Use DionysosLTL's generic simulator signature:
    #   simulate_buchi_policy(ts, aut, labels, policy, next_q; x0, q0, max_steps, stop_at_accept)
    _t_sim = @elapsed begin
        global xs, iqs, us, status = simulate_buchi_policy(
            ts, aut, labels_ts, polP, next_qP;
            x0 = x0_abs,
            q0 = q0_aut,
            max_steps = 1200,
            stop_at_accept = true,
        )
    end
    println("Büchi policy simulation time: ", _t_sim, " seconds")
    println("Total (labels + synth + sim) time: ", (_t_labels + _t_synth + _t_sim), " seconds")

    println("Büchi simulation status: ", status)
    println("inputs along Büchi run: ", us)

    # Plot abstract path in 2D using cell centers
    x1_path = Float64[]
    x2_path = Float64[]
    for s in xs
        x1c, x2c = center2d(s)
        push!(x1_path, x1c)
        push!(x2_path, x2c)
    end

    fig = plot_environment()
    xlims!(fig, (0.0, 10.0))
    ylims!(fig, (0.0, 10.0)) 
    plot!(fig, x1_path, x2_path;
          linewidth = 2,
          color = :red,
          marker = :circle,
          ms = 3,
          label = "Safe Trajectory")
    scatter!(fig, [x1_path[1]], [x2_path[1]]; marker = :diamond, ms = 9, color = :red, label = "start")
    scatter!(fig, [x1_path[end]], [x2_path[end]]; marker = :star5, ms = 9, color = :red, label = "end" ) # "end")
    display(fig)
    savefig(fig, "BNLX2026.png")
end # let
elseif RUN_STATIC_DEMO
let

    # --------------------------------------------------------------------------
    # TS-only reach-avoid synthesis (robust semantics)
    # Goal: reach :green while avoiding :obs
    # Robust semantics means: choose u such that ALL successors stay in the winning set.
    # --------------------------------------------------------------------------

    nx_ts = L2A.nstates(ts)
    nu_ts = L2A.ninputs(ts)

    # Goal/Avoid sets
    Goal = falses(nx_ts)
    Avoid = falses(nx_ts)
    for x in 1:nx_ts
        lbl = labels_ts[x]
        Goal[x]  = (:green in lbl)
        Avoid[x] = (:obs   in lbl)
    end
    println("TS-only: goal(:green)=", count(Goal), ", avoid(:obs)=", count(Avoid), ", total=", nx_ts)

    # Controllable predecessor on TS (robust): ∃u : post(x,u) ⊆ S
    function pre_ts_robust(ts::TSFromAdj, S::BitVector)
        nx = L2A.nstates(ts)
        nu = L2A.ninputs(ts)
        Pre = falses(nx)
        @inbounds for x in 1:nx
            ok_u = false
            for u in 1:nu
                succs = L2A.post(ts, x, u)
                isempty(succs) && continue
                all_in = true
                for xp in succs
                    if !S[xp]
                        all_in = false
                        break
                    end
                end
                if all_in
                    ok_u = true
                    break
                end
            end
            Pre[x] = ok_u
        end
        return Pre
    end

    function compute_reachavoid_robust(ts::TSFromAdj, Goal::BitVector, Avoid::BitVector)
        nx = L2A.nstates(ts)

        # Winning set initialization: goal but not avoid
        Win = falses(nx)
        @inbounds for x in 1:nx
            Win[x] = Goal[x] && !Avoid[x]
        end

        # Layer-to-goal (used for policy extraction)
        layer = fill(typemax(Int), nx)
        @inbounds for x in 1:nx
            if Win[x]
                layer[x] = 0
            end
        end

        # Backward fixed-point: Win = μZ. (Goal \ Avoid) ∪ Pre(Z) \ Avoid
        changed = true
        k = 0
        while changed
            changed = false
            k += 1
            PreW = pre_ts_robust(ts, Win)
            @inbounds for x in 1:nx
                if !Win[x] && PreW[x] && !Avoid[x]
                    Win[x] = true
                    layer[x] = k
                    changed = true
                end
            end
        end

        return Win, layer
    end

    Win_ts, layer_ts = compute_reachavoid_robust(ts, Goal, Avoid)

    println("TS-only: winning states (robust)=", count(Win_ts), " / ", nx_ts)

    x0_forced = ts.X0[1]
    println("Forced start x0_abs_start=", x0_forced,
            " is winning? robust=", Win_ts[x0_forced])

    # Extract a memoryless robust policy on TS
    policy_ts = fill(1, nx_ts)
    @inbounds for x in 1:nx_ts
        Win_ts[x] || continue
        if Goal[x]
            policy_ts[x] = 1
            continue
        end

        best_u = 1
        best_score = typemax(Int)
        best_tie = (typemax(Int), typemax(Float64), typemax(Float64), typemax(Int))

        for u in 1:nu_ts
            succs = L2A.post(ts, x, u)
            isempty(succs) && continue

            ok = true
            score = 0
            for xp in succs
                if !Win_ts[xp]
                    ok = false
                    break
                end
                score = max(score, layer_ts[xp])
            end
            ok || continue

            u_vec = input_from_index(u, abstract_system)
            has_move = any(xp -> xp != x, succs)
            tie = (u_vec[1] > 0 ? 0 : 1, abs(u_vec[2]), -abs(u_vec[1]), has_move ? 0 : 1)

            if (score < best_score) || (score == best_score && tie < best_tie)
                best_score = score
                best_tie = tie
                best_u = u
            end
        end

        policy_ts[x] = best_u
    end

    # Pick a winning initial abstract state from X0
    X0 = L2A.initset(ts)
    isempty(X0) && error("initset(ts) is empty; increase eps in X_start or change inclusion mode when building X0")

    idx0 = findfirst(s0 -> Win_ts[s0], X0)

    if idx0 === nothing
        @warn "Forced initial TS state is not winning under robust semantics; plotting diagnostics and stopping." length_X0=length(X0)

        n_obs0   = count(s0 -> (:obs in labels_ts[s0]), X0)
        n_green0 = count(s0 -> (:green in labels_ts[s0]), X0)
        println("Initial-set label counts: obs=", n_obs0, ", green=", n_green0, ", total=", length(X0))

        # Nearest robust-winning state distance in xy
        x0_xy = SVector(x0_cont[1], x0_cont[2])
        best_s = 0
        best_d = Inf
        for s in 1:nx_ts
            Win_ts[s] || continue
            x1c, x2c = center2d(s)
            d = norm(SVector(x1c, x2c) - x0_xy)
            if d < best_d
                best_d = d
                best_s = s
            end
        end
        if best_s == 0
            println("No robust-winning states exist at all (Win_ts empty).")
        else
            println("Nearest robust-winning state to x0_cont: s=", best_s, ", xy-distance=", best_d)
        end

        # Plot a local neighborhood of winning states around the forced start
        fig_dbg = plot_environment()
        x1s, x2s = center2d(x0_forced)
        scatter!(fig_dbg, [x1s], [x2s]; marker=:star5, ms=10, color=:black, label="forced start cell")
        scatter!(fig_dbg, [x0_cont[1]], [x0_cont[2]]; marker=:circle, ms=6, color=:gray, label="x0_cont")

        r = 3.0 * hx[1]
        wx = Float64[]; wy = Float64[]
        for s in 1:nx_ts
            Win_ts[s] || continue
            x1c, x2c = center2d(s)
            if norm(SVector(x1c, x2c) - SVector(x1s, x2s)) <= r
                push!(wx, x1c); push!(wy, x2c)
            end
        end
        if !isempty(wx)
            scatter!(fig_dbg, wx, wy; ms=3, color=:magenta, label= false ) # "robust-winning nearby")
        end
        display(fig_dbg)

        println("Stopping: the forced start cell is losing under robust reach-avoid semantics.")
        println("To make it robust-winning, refine the abstraction (smaller hx / tighter growth bound / possibly smaller DT).")

    else
        x0_abs = X0[idx0]
        println("Chosen winning initial abstract state (TS-only): ", x0_abs, ", label=", labels_ts[x0_abs])

        # Simulate a representative TS run
        ts_path_off = Int[x0_abs]
        u_path_off  = Int[]

        max_steps = 400
        status = :max_steps
        for step in 1:max_steps
            x = ts_path_off[end]
            Goal[x] && (status = :hit_goal; break)
            u = policy_ts[x]
            push!(u_path_off, u)
            succs = L2A.post(ts, x, u)
            isempty(succs) && (status = :deadend; break)

            # representative successor: prefer movement and smallest layer
            best = nothing
            best_layer = typemax(Int)
            for cand in succs
                Avoid[cand] && continue
                Win_ts[cand] || continue
                cand_layer = layer_ts[cand]
                if cand != x && cand_layer < best_layer
                    best = cand
                    best_layer = cand_layer
                end
            end
            if best === nothing
                for cand in succs
                    Avoid[cand] && continue
                    Win_ts[cand] || continue
                    cand_layer = layer_ts[cand]
                    if cand_layer < best_layer
                        best = cand
                        best_layer = cand_layer
                    end
                end
            end
            best === nothing && (status = :left_winning_set; break)

            push!(ts_path_off, best)
        end

        println("TS-only simulation status: ", status)
        println("inputs along TS-only run: ", u_path_off)

        # Build abstract path in (x1,x2)
        x1_path = Float64[]
        x2_path = Float64[]
        for s in ts_path_off
            x1c, x2c = center2d(s)
            push!(x1_path, x1c)
            push!(x2_path, x2c)
        end

        fig = plot_environment()
        plot!(fig, x1_path, x2_path;
              linewidth = 3,
              color = :red,
              marker = :circle,
              ms = 4,
              label = "TS-only abstract path")
        scatter!(fig, [x1_path[1]], [x2_path[1]]; marker = :star5, ms = 9, color = :red, label = false ) # "start")
        scatter!(fig, [x1_path[end]], [x2_path[end]]; marker = :diamond, ms = 9, color = :red, label= false ) # "end")
        display(fig)

        println("Plotted TS-only abstract path (", length(x1_path), " points).")
    end

end # let
end # if RUN_STATIC_DEMO / RUN_BUCHI_GAME

# --------------------------------------------------------------------------
# Genie web frontend.
# --------------------------------------------------------------------------

include("Genie_interface.jl")

const RUN_GENIE_FRONTEND = true
if RUN_GENIE_FRONTEND
    global GENIE_SERVER = start_genie_interface(; async = true)
end

