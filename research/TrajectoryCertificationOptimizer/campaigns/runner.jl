# Shared campaign runner for the robustness / ablation campaigns of plan.md §6b.
#
# A campaign is a grid of configurations, each evaluated over many rng seeds — the
# pipeline is stochastic (MPPI), so a configuration's result is a success rate and
# quantiles, never a single run. Campaign scripts stay ~50 lines: they define the
# configs and a `run_one(config, rng) -> NamedTuple` and call `run_campaign`.
#
# Output: one CSV row per (config, seed) under `campaigns/results/<name>.csv`
# (the raw material), plus a printed per-config aggregate table (the summary).
# Plots and animations are regenerated, never committed.
#
# Run with the bench environment (has Dionysos, Clarabel, Symbolics, CSV):
#     julia --project=bench research/TrajectoryCertificationOptimizer/campaigns/<campaign>.jl

module CampaignRunner

using Printf
using Statistics
using Random
import CSV
import DataFrames

export run_campaign

# NaN entries are legitimate (a failed run has no volume/time to report) —
# quantiles are taken over the defined values only.
function _quantiles(v)
    defined = filter(!isnan, v)
    isempty(defined) && return (NaN, NaN, NaN)
    return (quantile(defined, 0.25), median(defined), quantile(defined, 0.75))
end

"""
    run_campaign(; name, configs, run_one, nseeds = 20, results_dir = nothing)

Evaluate every `config ∈ configs` over `nseeds` seeded runs of
`run_one(config, rng)::NamedTuple`. The returned tuple must contain
`success::Bool`; every other `Real` field is aggregated (median and interquartile
range over the seeds). Each config needs a `label::String` field.

Writes the per-run rows to `<results_dir>/<name>.csv` and prints the aggregate
table. Returns the per-run rows.
"""
function run_campaign(;
    name::AbstractString,
    configs,
    run_one,
    nseeds::Integer = 20,
    results_dir = nothing,
)
    rows = NamedTuple[]
    aggregates = NamedTuple[]

    for config in configs
        outs = NamedTuple[]
        for s in 1:nseeds
            rng = Random.MersenneTwister(1000 + s)
            out = run_one(config, rng)
            push!(outs, out)
            push!(rows, (; label = config.label, seed = s, out...))
        end

        numeric = [k for k in keys(first(outs)) if k != :success && first(outs)[k] isa Real]
        agg = Dict{Symbol, Any}(
            :label => config.label,
            :success_rate => mean(o.success for o in outs),
        )
        for k in numeric
            q25, med, q75 = _quantiles([Float64(o[k]) for o in outs])
            agg[Symbol(k, :_median)] = med
            agg[Symbol(k, :_iqr)] = q75 - q25
        end
        push!(aggregates, NamedTuple(agg))
    end

    _print_aggregates(name, aggregates, nseeds)

    dir = results_dir === nothing ? joinpath(@__DIR__, "results") : results_dir
    isdir(dir) || mkpath(dir)
    path = joinpath(dir, "$name.csv")
    CSV.write(path, DataFrames.DataFrame(rows))
    println("Wrote $path")

    return rows
end

function _print_aggregates(name, aggregates, nseeds)
    println("\n" * repeat("=", 84))
    println("CAMPAIGN $name  ($nseeds seeds per config)")
    println(repeat("=", 84))
    for agg in aggregates
        @printf("%-24s success %5.1f%%", agg.label, 100 * agg.success_rate)
        for k in sort(collect(keys(agg)))
            (k == :label || k == :success_rate) && continue
            endswith(String(k), "_iqr") && continue
            base = replace(String(k), "_median" => "")
            iqr = get(agg, Symbol(base, :_iqr), NaN)
            @printf("   %s %.3g (iqr %.2g)", base, agg[k], iqr)
        end
        println()
    end
    println(repeat("=", 84))
    return nothing
end

end # module
