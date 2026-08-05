#!/usr/bin/env julia

"""
Dionysos cross-solver benchmark suite
=====================================
Runs each (solver family, problem) case from `definitions.jl`, timing the full
setup → abstraction build → controller synthesis, and reports a comparison table
plus a CSV. Every case is wrapped in try/catch, so one failure never aborts the run.

Usage:
    julia --project=bench bench/solver_suite.jl [trials] [filter]

Arguments:
    trials  : runs per case; the fastest run is reported (default: 1).
              Use >= 2 to discard Julia's one-time compilation latency.
    filter  : only run cases whose "method / problem" label contains this substring
              (case-insensitive), e.g. `Gol` or `UniformGrid` (default: run all).

Output:
    bench/results/solver_suite.csv
"""

using Printf
using Statistics
import MathOptInterface as MOI

include(joinpath(@__DIR__, "definitions.jl"))

import CSV
import DataFrames

function _termination_status(optimizer)
    # MIQP families expose an MOI status; abstraction solvers don't, so a completed
    # `optimize!` reports "completed".
    return try
        string(MOI.get(optimizer, MOI.TerminationStatus()))
    catch
        "completed"
    end
end

"""Run one case `trials` times; return the fastest timing plus status/allocations."""
function run_case(case; trials::Int = 1)
    best_time = Inf
    best_bytes = 0
    status = "—"
    for _ in 1:trials
        result = @timed case.run()
        if result.time < best_time
            best_time = result.time
            best_bytes = result.bytes
            status = _termination_status(result.value)
        end
    end
    return (;
        ok = true,
        status = status,
        time_s = best_time,
        alloc_bytes = best_bytes,
        error = "",
    )
end

function format_bytes(bytes::Real)
    units = ("B", "KB", "MB", "GB", "TB")
    value = float(bytes)
    i = 1
    while value >= 1024 && i < length(units)
        value /= 1024
        i += 1
    end
    return @sprintf("%.2f %s", value, units[i])
end

function print_table(rows)
    println("\n" * repeat("=", 104))
    println("CROSS-SOLVER BENCHMARK SUITE")
    println(repeat("=", 104))
    @printf(
        "%-28s %-26s %-22s %11s %13s\n",
        "Method",
        "Problem",
        "Status",
        "Time (s)",
        "Allocations"
    )
    println(repeat("-", 104))
    for r in rows
        timestr = isnan(r.time_s) ? "—" : @sprintf("%.3f", r.time_s)
        @printf(
            "%-28s %-26s %-22s %11s %13s\n",
            r.method,
            r.problem,
            r.ok ? r.status : "ERROR",
            timestr,
            r.ok ? format_bytes(r.alloc_bytes) : "—",
        )
    end
    println(repeat("=", 104))
    failures = count(r -> !r.ok, rows)
    failures == 0 || println("$failures case(s) failed — see the error column in the CSV.")
    return nothing
end

function write_csv(rows)
    dir = joinpath(@__DIR__, "results")
    isdir(dir) || mkpath(dir)
    path = joinpath(dir, "solver_suite.csv")
    df = DataFrames.DataFrame(;
        method = [r.method for r in rows],
        problem = [r.problem for r in rows],
        status = [r.ok ? r.status : "ERROR" for r in rows],
        ok = [r.ok for r in rows],
        time_s = [r.time_s for r in rows],
        alloc_bytes = [r.alloc_bytes for r in rows],
        error = [r.error for r in rows],
    )
    CSV.write(path, df)
    println("Wrote $path")
    return path
end

function main()
    trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    filter = length(ARGS) >= 2 ? lowercase(ARGS[2]) : ""

    cases = [
        c for c in CASES if
        isempty(filter) || occursin(filter, lowercase("$(c.method) / $(c.problem)"))
    ]
    isempty(cases) && (println("No cases match filter \"$filter\"."); return nothing)

    println("Running $(length(cases)) case(s), $trials trial(s) each.")
    rows = NamedTuple[]
    for case in cases
        @printf("→ %-28s %-26s ... ", case.method, case.problem)
        flush(stdout)
        row = try
            r = run_case(case; trials = trials)
            @printf("%s  %.3fs\n", r.status, r.time_s)
            (;
                method = case.method,
                problem = case.problem,
                ok = true,
                status = r.status,
                time_s = r.time_s,
                alloc_bytes = r.alloc_bytes,
                error = "",
            )
        catch e
            msg = sprint(showerror, e)
            println("FAILED: ", first(split(msg, '\n')))
            (;
                method = case.method,
                problem = case.problem,
                ok = false,
                status = "ERROR",
                time_s = NaN,
                alloc_bytes = 0,
                error = replace(msg, '\n' => " | "),
            )
        end
        push!(rows, row)
    end

    print_table(rows)
    write_csv(rows)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
