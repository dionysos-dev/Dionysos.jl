#!/usr/bin/env julia
"""
Assemble benchmark results from JSON metrics files and produce a LaTeX report.

Usage:
    julia --project=<PROJECT> generate_report.jl [results_dir] [report_dir]

Reads all `metrics_*.json` files under results_dir (recursively) and generates
a LaTeX comparison report in report_dir.
"""

using JSON
using Dates
using Printf

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function find_json_files(dir::String)
    files = String[]
    for (root, _, filenames) in walkdir(dir)
        for f in filenames
            if startswith(f, "metrics_") && endswith(f, ".json")
                push!(files, joinpath(root, f))
            end
        end
    end
    return files
end

function load_results(results_dir::String)
    results = Dict{String, Vector{Dict{String, Any}}}()
    for f in find_json_files(results_dir)
        data = JSON.parsefile(f)
        mode = get(data, "mode", "unknown")
        if !haskey(results, mode)
            results[mode] = []
        end
        push!(results[mode], data)
    end
    return results
end

function latest_result(runs::Vector{Dict{String, Any}})
    sort!(runs; by = r -> get(r, "timestamp", ""), rev = true)
    return runs[1]
end

safe_get(d::Dict, key::String, default = "—") = get(d, key, default)

function fmt_time(sec)
    return sec isa Number ? @sprintf("%.4f", sec) : string(sec)
end

function fmt_float(val; digits = 2)
    return val isa Number ? @sprintf("%.*f", digits, val) : string(val)
end

function escape_latex(s::String)
    s = replace(s, "\\" => "\\textbackslash{}")
    s = replace(s, "&" => "\\&")
    s = replace(s, "%" => "\\%")
    s = replace(s, "\$" => "\\\$")
    s = replace(s, "#" => "\\#")
    s = replace(s, "_" => "\\_")
    s = replace(s, "{" => "\\{")
    s = replace(s, "}" => "\\}")
    s = replace(s, "~" => "\\textasciitilde{}")
    s = replace(s, "^" => "\\textasciicircum{}")
    return s
end

# ---------------------------------------------------------------------------
# Log parsing — extract structured info from run.log files
# ---------------------------------------------------------------------------

"""
    parse_log_file(log_path) -> Dict{String, Any}

Parse a `run.log` file and extract known key-value patterns.
Returns a dictionary with parsed values (keys that are not found are omitted).
"""
function parse_log_file(log_path::String)
    isfile(log_path) || return Dict{String, Any}()
    text = read(log_path, String)
    parsed = Dict{String, Any}()

    # Patterns from example output
    patterns = [
        (
            r"Time to construct the abstraction:\s*([\d.eE+-]+)",
            "abstraction_time_sec",
            Float64,
        ),
        (r"Time to solve the abstract problem:\s*([\d.eE+-]+)", "solve_time_sec", Float64),
        (r"n_state:\s*(\d+)", "n_state", Int),
        (r"n_input:\s*(\d+)", "n_input", Int),
        (r"n_pos:\s*(\d+)", "n_pos", Int),
        (r"n_vel:\s*(\d+)", "n_vel", Int),
    ]

    # TIMING patterns from instrumented code
    timing_patterns = [
        (r"TIMING loading_time:\s*([\d.eE+-]+)", "loading_time_sec", Float64),
        (r"TIMING system_setup_time:\s*([\d.eE+-]+)", "system_setup_time_sec", Float64),
        (
            r"TIMING pre_abstraction_time:\s*([\d.eE+-]+)",
            "pre_abstraction_time_sec",
            Float64,
        ),
        (
            r"TIMING post_abstraction_time:\s*([\d.eE+-]+)",
            "post_abstraction_time_sec",
            Float64,
        ),
        (r"TIMING rbd_simulate_time:\s*([\d.eE+-]+)", "rbd_simulate_time_sec", Float64),
        (r"TIMING rbd_call_count:\s*(\d+)", "rbd_call_count", Int),
    ]

    # Patterns from wrapper header
    header_patterns = [
        (r"NPARTS:\s*(\d+)", "log_nparts", Int),
        (r"Threads:\s*(\d+)", "log_threads", Int),
        (r"PID:\s*(\d+)", "log_pid", Int),
        (r"Host:\s*(\S+)", "log_host", String),
        (r"Mode:\s*(\S+)", "log_mode", String),
    ]

    # Patterns from wrapper footer
    footer_patterns = [
        (r"Wall-clock:\s*([\d.]+)\s*s", "log_wall_clock_sec", Float64),
        (r"@timed:\s*([\d.]+)\s*s", "log_timed_sec", Float64),
        (r"Alloc:\s*([\d.]+)\s*MB", "log_alloc_MB", Float64),
        (r"GC time:\s*([\d.]+)\s*s", "log_gc_time_sec", Float64),
    ]

    for (regex, key, T) in vcat(patterns, timing_patterns, header_patterns, footer_patterns)
        m = match(regex, text)
        if m !== nothing
            try
                parsed[key] = T == String ? m[1] : parse(T, m[1])
            catch
                parsed[key] = m[1]
            end
        end
    end

    # Count "addprocs" evidence: look for worker-related messages
    if occursin(r"Workers removed", text)
        parsed["had_workers"] = true
    end

    # Parse per-worker timing lines
    worker_times = Dict{String, Any}[]
    for m in eachmatch(
        r"TIMING worker_time worker_id=(\d+) partition=(\d+) elapsed=([\d.eE+-]+) n_states=(\d+) n_transitions=(\d+)",
        text,
    )
        push!(
            worker_times,
            Dict{String, Any}(
                "worker_id" => parse(Int, m[1]),
                "partition" => parse(Int, m[2]),
                "elapsed" => parse(Float64, m[3]),
                "n_states" => parse(Int, m[4]),
                "n_transitions" => parse(Int, m[5]),
            ),
        )
    end
    if !isempty(worker_times)
        parsed["worker_times"] = worker_times
    end

    # Parse per-worker RBD timing lines
    rbd_worker_times = Dict{String, Any}[]
    for m in
        eachmatch(r"TIMING rbd_worker worker_id=(\d+) time=([\d.eE+-]+) calls=(\d+)", text)
        push!(
            rbd_worker_times,
            Dict{String, Any}(
                "worker_id" => parse(Int, m[1]),
                "time_sec" => parse(Float64, m[2]),
                "call_count" => parse(Int, m[3]),
            ),
        )
    end
    if !isempty(rbd_worker_times)
        parsed["rbd_worker_times"] = rbd_worker_times
    end

    # Parse per-thread timing lines
    thread_times = Dict{String, Any}[]
    for m in eachmatch(
        r"TIMING thread_time thread_id=(\d+) elapsed=([\d.eE+-]+) n_transitions=(\d+)",
        text,
    )
        push!(
            thread_times,
            Dict{String, Any}(
                "thread_id" => parse(Int, m[1]),
                "elapsed" => parse(Float64, m[2]),
                "n_transitions" => parse(Int, m[3]),
            ),
        )
    end
    if !isempty(thread_times)
        parsed["thread_times"] = thread_times
    end

    return parsed
end

"""
    load_log_data(results_dir) -> Dict{String, Dict{String, Any}}

Find and parse run.log files for each mode subdirectory.
"""
function load_log_data(results_dir::String)
    logs = Dict{String, Dict{String, Any}}()
    for mode in ["serial", "threaded", "distributed", "hybrid"]
        log_path = joinpath(results_dir, mode, "run.log")
        if isfile(log_path)
            logs[mode] = parse_log_file(log_path)
        end
    end
    # Also scan scaling subdirectories (distributed_n*, threaded_n*)
    if isdir(results_dir)
        for entry in readdir(results_dir)
            m = match(r"^(distributed|threaded)_n(\d+)$", entry)
            m === nothing && continue
            log_path = joinpath(results_dir, entry, "run.log")
            if isfile(log_path)
                logs[entry] = parse_log_file(log_path)
            end
        end
    end
    return logs
end

# ---------------------------------------------------------------------------
# LaTeX generation
# ---------------------------------------------------------------------------

function generate_latex(
    results::Dict{String, Vector{Dict{String, Any}}},
    report_dir::String;
    log_data::Dict{String, Dict{String, Any}} = Dict{String, Dict{String, Any}}(),
)
    modes = ["serial", "threaded", "distributed", "hybrid"]
    data = Dict{String, Dict{String, Any}}()

    for m in modes
        if haskey(results, m) && !isempty(results[m])
            data[m] = latest_result(results[m])
        end
    end

    if isempty(data)
        println("ERROR: No results found. Cannot generate report.")
        return nothing
    end

    tex = IOBuffer()

    # ── Preamble ──
    example_name =
        escape_latex(string(safe_get(first(values(data)), "example", "Dionysos")))

    println(
        tex,
        raw"""
\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2.5cm]{geometry}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage{colortbl}
\usepackage{hyperref}
\usepackage{graphicx}
\usepackage{float}
\usepackage{caption}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{enumitem}
\usepackage{fancyhdr}

\setlength{\headheight}{14pt}
\addtolength{\topmargin}{-2pt}

\definecolor{headerblue}{RGB}{41,65,122}
\definecolor{lightgray}{RGB}{240,240,240}
\definecolor{successgreen}{RGB}{34,139,34}
\definecolor{warnorange}{RGB}{200,120,0}

\pagestyle{fancy}
\fancyhf{}
\rhead{Dionysos Parallel Benchmarks}
\lhead{\today}
\cfoot{\thepage}

\hypersetup{
    colorlinks=true,
    linkcolor=headerblue,
    urlcolor=headerblue,
}
""",
    )
    println(
        tex,
        "\\title{\\textcolor{headerblue}{\\textbf{Dionysos Parallel Benchmark Report}}\\\\[0.5em]",
    )
    println(
        tex,
        "       \\large $(example_name): Serial vs.\\ Threaded vs.\\ Distributed vs.\\ Hybrid}",
    )
    println(tex, "\\author{Auto-generated by \\texttt{generate\\_report.jl}}")
    println(tex, "\\date{", escape_latex(string(Dates.now())), "}")
    println(
        tex,
        raw"""
\begin{document}
\maketitle
\tableofcontents
\newpage
""",
    )

    # ── Section 1: Executive Summary ──
    println(
        tex,
        raw"""
\section{Executive Summary}
This report compares four parallelism strategies for Dionysos abstraction-based controller synthesis:
\begin{enumerate}[leftmargin=2em]
  \item \textbf{Serial} — single thread, single process (baseline).
  \item \textbf{Multithreaded} — multiple Julia threads within one process (\texttt{Threads.@threads}).
  \item \textbf{Distributed} — multiple Julia processes (workers) via \texttt{Distributed.jl} + \texttt{pmap}.
  \item \textbf{Hybrid} — multiple workers, each using multiple threads.
\end{enumerate}
The abstraction construction (computing the symbolic transition system) is the parallelized phase.
Each mode is controlled via the \texttt{DIONYSOS\_DISTRIBUTED} and \texttt{DIONYSOS\_THREADED} environment variables.
""",
    )

    # ── Section 2: Environment ──
    println(tex, raw"\section{Environment}")
    sample = first(values(data))
    println(tex, raw"\begin{table}[H]")
    println(tex, raw"\centering")
    println(tex, raw"\caption{Runtime Environment}")
    println(tex, raw"\begin{tabular}{ll}")
    println(tex, raw"\toprule")
    println(tex, raw"\textbf{Property} & \textbf{Value} \\\\")
    println(tex, raw"\midrule")
    println(
        tex,
        "Julia version & \\texttt{",
        escape_latex(string(safe_get(sample, "julia_version"))),
        "} \\\\",
    )
    println(
        tex,
        "Hostname & \\texttt{",
        escape_latex(string(safe_get(sample, "hostname"))),
        "} \\\\",
    )
    println(tex, "Total memory & ", safe_get(sample, "total_memory_GB"), " GB \\\\")
    println(
        tex,
        "Example & \\texttt{",
        escape_latex(string(safe_get(sample, "example"))),
        "} \\\\",
    )

    for var in [
        "SLURM_JOB_ID",
        "SLURM_NNODES",
        "SLURM_NTASKS",
        "SLURM_CPUS_PER_TASK",
        "SLURM_NODELIST",
    ]
        val = safe_get(sample, var, nothing)
        if val !== nothing
            println(
                tex,
                escape_latex(var),
                " & \\texttt{",
                escape_latex(string(val)),
                "} \\\\",
            )
        end
    end
    println(tex, raw"\bottomrule")
    println(tex, raw"\end{tabular}")
    println(tex, raw"\end{table}")

    # ── Section 3: Problem Characteristics (from logs) ──
    if !isempty(log_data)
        # Use any available log to get problem dimensions
        sample_log = first(values(log_data))
        has_dims =
            any(haskey(sample_log, k) for k in ["n_state", "n_input", "n_pos", "n_vel"])
        if has_dims
            println(tex, raw"\section{Problem Characteristics}")
            println(tex, raw"\begin{table}[H]")
            println(tex, raw"\centering")
            println(tex, raw"\caption{Problem dimensions (from example output)}")
            println(tex, raw"\begin{tabular}{ll}")
            println(tex, raw"\toprule")
            println(tex, raw"\textbf{Property} & \textbf{Value} \\\\")
            println(tex, raw"\midrule")
            for (key, label) in [
                ("n_state", "State dimension"),
                ("n_input", "Input dimension"),
                ("n_pos", "Position DoF"),
                ("n_vel", "Velocity DoF"),
            ]
                val = get(sample_log, key, nothing)
                val !== nothing && println(tex, "$label & \\texttt{$val} \\\\")
            end
            println(tex, raw"\bottomrule")
            println(tex, raw"\end{tabular}")
            println(tex, raw"\end{table}")
        end
    end

    # ── Section 4: Combined Results & Abstraction Table ──
    println(
        tex,
        raw"""
\section{Performance Comparison}
\subsection{Combined Summary}
""",
    )

    serial_time = nothing
    if haskey(data, "serial")
        serial_time = get(data["serial"], "elapsed_sec", nothing)
    end

    serial_abs_time = nothing
    abs_times = Dict{String, Float64}()
    for m in modes
        lg = get(log_data, m, Dict())
        abt = get(lg, "abstraction_time_sec", nothing)
        if abt !== nothing
            abs_times[m] = abt
            m == "serial" && (serial_abs_time = abt)
        end
    end

    println(tex, raw"\begin{table}[H]")
    println(tex, raw"\centering")
    println(tex, raw"\caption{Benchmark and abstraction results — all modes}")
    println(tex, raw"\begin{tabular}{l r r r r r r r r}")
    println(tex, raw"\toprule")
    println(
        tex,
        raw"\textbf{Mode} & \textbf{Total (s)} & \textbf{Speedup} & \textbf{Abstr.\ (s)} & \textbf{Abstr.\ Speedup} & \textbf{Threads} & \textbf{Workers} & \textbf{Parts} & \textbf{Alloc (MB)} \\\\",
    )
    println(tex, raw"\midrule")

    for m in modes
        haskey(data, m) || continue
        d = data[m]
        t = get(d, "elapsed_sec", 0.0)
        threads = get(d, "julia_threads", 1)
        nprocs_actual = get(
            d,
            "example_nprocs",
            get(get(log_data, m, Dict()), "log_nparts", get(d, "dionysos_nparts", 1)),
        )
        parts = get(d, "dionysos_nparts", 1)
        alloc = get(d, "alloc_MB", "—")

        speedup = "1.00"
        if serial_time !== nothing && serial_time > 0 && t > 0
            speedup = fmt_float(serial_time / t)
        end

        abt_str = "—"
        abs_sp_str = "—"
        if haskey(abs_times, m)
            abt_str = fmt_time(abs_times[m])
            if serial_abs_time !== nothing && serial_abs_time > 0 && abs_times[m] > 0
                abs_sp_str = fmt_float(serial_abs_time / abs_times[m])
            end
        end

        mode_label = Dict(
            "serial" => "Serial",
            "threaded" => "Threaded",
            "distributed" => "Distributed",
            "hybrid" => "Hybrid",
        )[m]
        println(
            tex,
            "$mode_label & $(fmt_time(t)) & $(speedup)\\texttimes & $abt_str & $(abs_sp_str)\\texttimes & $threads & $nprocs_actual & $parts & $alloc \\\\",
        )
    end

    println(tex, raw"\bottomrule")
    println(tex, raw"\end{tabular}")
    println(tex, raw"\end{table}")

    # ── Timing Breakdown Table ──
    has_breakdown = any(haskey(get(log_data, m, Dict()), "loading_time_sec") for m in modes)
    if has_breakdown
        println(tex, raw"\subsection{Timing Breakdown}")
        println(tex, raw"\begin{table}[H]")
        println(tex, raw"\centering")
        println(tex, raw"\caption{Time breakdown per phase (seconds)}")

        available_modes = [m for m in modes if haskey(data, m)]
        ncols = length(available_modes) + 1
        col_spec = "l" * repeat(" r", length(available_modes))
        println(tex, "\\begin{tabular}{$col_spec}")
        println(tex, raw"\toprule")
        mode_headers = join(
            [
                "\\textbf{$(Dict("serial"=>"Serial","threaded"=>"Threaded","distributed"=>"Distributed","hybrid"=>"Hybrid")[m])}"
                for m in available_modes
            ],
            " & ",
        )
        println(tex, "\\textbf{Phase} & $mode_headers \\\\")
        println(tex, raw"\midrule")

        breakdown_keys = [
            ("loading_time_sec", "Loading \\& imports"),
            ("system_setup_time_sec", "System setup (RBD init)"),
            ("pre_abstraction_time_sec", "Pre-abstraction total"),
            ("abstraction_time_sec", "Abstraction construction"),
            ("rbd_simulate_time_sec", "RigidBodyDynamics (cumul.)"),
            ("rbd_call_count", "RBD call count"),
            ("post_abstraction_time_sec", "Post-abstraction total"),
            ("log_wall_clock_sec", "Wall-clock total"),
        ]

        for (key, label) in breakdown_keys
            vals = String[]
            for m in available_modes
                lg = get(log_data, m, Dict())
                val = get(lg, key, nothing)
                if val === nothing
                    push!(vals, "—")
                elseif val isa AbstractFloat
                    push!(vals, fmt_time(val))
                else
                    push!(vals, string(val))
                end
            end
            println(tex, "$label & $(join(vals, " & ")) \\\\")
        end

        println(tex, raw"\bottomrule")
        println(tex, raw"\end{tabular}")
        println(tex, raw"\end{table}")
    end

    # ── Per-Worker / Per-Thread Timing ──
    has_worker_or_thread_times = any(
        haskey(get(log_data, m, Dict()), "worker_times") ||
        haskey(get(log_data, m, Dict()), "thread_times") ||
        haskey(get(log_data, m, Dict()), "rbd_worker_times") for m in modes
    )
    if has_worker_or_thread_times
        println(tex, raw"\subsection{Per-Worker / Per-Thread Timing}")

        for m in modes
            lg = get(log_data, m, Dict())
            wt = get(lg, "worker_times", nothing)
            tt = get(lg, "thread_times", nothing)

            mode_label = Dict(
                "serial" => "Serial",
                "threaded" => "Threaded",
                "distributed" => "Distributed",
                "hybrid" => "Hybrid",
            )[m]

            if wt !== nothing && !isempty(wt)
                println(tex, "\n\\subsubsection{$mode_label — Worker Times}")
                println(tex, raw"\begin{table}[H]")
                println(tex, raw"\centering")
                println(tex, raw"\begin{tabular}{r r r r r}")
                println(tex, raw"\toprule")
                println(
                    tex,
                    raw"\textbf{Worker ID} & \textbf{Partition} & \textbf{Elapsed (s)} & \textbf{States} & \textbf{Transitions} \\\\",
                )
                println(tex, raw"\midrule")
                for w in wt
                    println(
                        tex,
                        "$(w["worker_id"]) & $(w["partition"]) & $(fmt_time(w["elapsed"])) & $(w["n_states"]) & $(w["n_transitions"]) \\\\",
                    )
                end
                println(tex, raw"\bottomrule")
                println(tex, raw"\end{tabular}")
                println(tex, raw"\end{table}")
            end

            rbd_wt = get(lg, "rbd_worker_times", nothing)
            if rbd_wt !== nothing && !isempty(rbd_wt)
                println(tex, "\n\\subsubsection{$mode_label — RBD Per-Worker Times}")
                println(tex, raw"\begin{table}[H]")
                println(tex, raw"\centering")
                println(tex, raw"\begin{tabular}{r r r}")
                println(tex, raw"\toprule")
                println(
                    tex,
                    raw"\textbf{Worker ID} & \textbf{RBD Time (s)} & \textbf{RBD Calls} \\\\",
                )
                println(tex, raw"\midrule")
                for w in rbd_wt
                    println(
                        tex,
                        "$(w["worker_id"]) & $(fmt_time(w["time_sec"])) & $(w["call_count"]) \\\\",
                    )
                end
                println(tex, raw"\bottomrule")
                println(tex, raw"\end{tabular}")
                println(tex, raw"\end{table}")
            end

            if tt !== nothing && !isempty(tt)
                println(tex, "\n\\subsubsection{$mode_label — Thread Times}")
                println(tex, raw"\begin{table}[H]")
                println(tex, raw"\centering")
                println(tex, raw"\begin{tabular}{r r r}")
                println(tex, raw"\toprule")
                println(
                    tex,
                    raw"\textbf{Thread ID} & \textbf{Elapsed (s)} & \textbf{Transitions} \\\\",
                )
                println(tex, raw"\midrule")
                for t in tt
                    println(
                        tex,
                        "$(t["thread_id"]) & $(fmt_time(t["elapsed"])) & $(t["n_transitions"]) \\\\",
                    )
                end
                println(tex, raw"\bottomrule")
                println(tex, raw"\end{tabular}")
                println(tex, raw"\end{table}")
            end
        end
    end

    # ── Section: Combined Detailed Metrics (4-column table) ──
    println(tex, raw"\subsection{Detailed Metrics (All Modes)}")

    available_modes = [m for m in modes if haskey(data, m)]
    if !isempty(available_modes)
        col_spec = "l" * repeat(" r", length(available_modes))
        mode_headers = join(
            [
                "\\textbf{$(Dict("serial"=>"Serial","threaded"=>"Threaded","distributed"=>"Distributed","hybrid"=>"Hybrid")[m])}"
                for m in available_modes
            ],
            " & ",
        )

        println(tex, raw"\begin{table}[H]")
        println(tex, raw"\centering")
        println(tex, raw"\caption{Detailed metrics comparison}")
        println(tex, "\\begin{tabular}{$col_spec}")
        println(tex, raw"\toprule")
        println(tex, "\\textbf{Metric} & $mode_headers \\\\")
        println(tex, raw"\midrule")

        combined_detail_keys = [
            ("elapsed_sec", "Wall-clock time (s)"),
            ("wall_clock_sec", "Total time incl.\\ setup (s)"),
            ("julia_threads", "Julia threads"),
            ("julia_nprocs", "Julia processes"),
            ("julia_nworkers", "Julia workers"),
            ("dionysos_nparts", "Partitions (NPARTS)"),
            ("alloc_MB", "Allocated memory (MB)"),
            ("gc_time_sec", "GC time (s)"),
            ("free_memory_GB", "Free memory (GB)"),
            ("total_memory_GB", "Total memory (GB)"),
        ]

        for (key, label) in combined_detail_keys
            vals = String[]
            for m in available_modes
                d = data[m]
                val = get(d, key, nothing)
                if val === nothing
                    push!(vals, "—")
                elseif val isa AbstractFloat
                    push!(vals, fmt_time(val))
                elseif val isa Bool
                    push!(vals, val ? "Yes" : "No")
                else
                    push!(vals, escape_latex(string(val)))
                end
            end
            println(tex, "$label & $(join(vals, " & ")) \\\\")
        end

        # Speedup row
        speed_vals = String[]
        for m in available_modes
            t = get(data[m], "elapsed_sec", 0.0)
            if serial_time !== nothing && serial_time > 0 && t > 0
                push!(speed_vals, "$(fmt_float(serial_time / t))\\texttimes")
            else
                push!(speed_vals, "—")
            end
        end
        println(tex, "Speedup vs.\\ serial & $(join(speed_vals, " & ")) \\\\")

        # Log-derived rows
        log_detail_keys = [
            ("abstraction_time_sec", "Abstraction time (s)"),
            ("loading_time_sec", "Loading time (s)"),
            ("system_setup_time_sec", "System setup time (s)"),
            ("pre_abstraction_time_sec", "Pre-abstraction time (s)"),
            ("rbd_simulate_time_sec", "RBD simulate time (s)"),
            ("rbd_call_count", "RBD call count"),
        ]

        println(tex, raw"\midrule")
        for (key, label) in log_detail_keys
            vals = String[]
            any_present = false
            for m in available_modes
                lg = get(log_data, m, Dict())
                val = get(lg, key, nothing)
                if val === nothing
                    push!(vals, "—")
                else
                    any_present = true
                    push!(vals, val isa AbstractFloat ? fmt_time(val) : string(val))
                end
            end
            any_present && println(tex, "$label & $(join(vals, " & ")) \\\\")
        end

        println(tex, raw"\bottomrule")
        println(tex, raw"\end{tabular}")
        println(tex, raw"\end{table}")
    end

    # ── Scaling Data & Plots ──
    scaling_dist = Dict{Int, Dict{String, Any}}()
    scaling_thr = Dict{Int, Dict{String, Any}}()
    for (key, lg) in log_data
        m = match(r"^distributed_n(\d+)$", key)
        if m !== nothing
            n = parse(Int, m[1])
            scaling_dist[n] = lg
        end
        m = match(r"^threaded_n(\d+)$", key)
        if m !== nothing
            n = parse(Int, m[1])
            scaling_thr[n] = lg
        end
    end

    has_scaling = !isempty(scaling_dist) || !isempty(scaling_thr)
    if has_scaling
        println(tex, raw"\subsection{Scaling Results}")

        # Generate scaling data files and PGF-free plots via LaTeX picture env
        # We write CSV files and include as tables; also generate gnuplot/pgfplots if available
        for (label, sdata, xlab) in
            [("Distributed", scaling_dist, "Workers"), ("Threaded", scaling_thr, "Threads")]
            isempty(sdata) && continue
            ns = sort(collect(keys(sdata)))

            println(tex, "\n\\subsubsection{$label Scaling}")
            println(tex, raw"\begin{table}[H]")
            println(tex, raw"\centering")
            println(tex, "\\caption{$label scaling data}")
            println(tex, raw"\begin{tabular}{r r r r r}")
            println(tex, raw"\toprule")
            println(
                tex,
                "\\textbf{$xlab} & \\textbf{Total (s)} & \\textbf{Abstr.\\ (s)} & \\textbf{Speedup (total)} & \\textbf{Speedup (abstr.)} \\\\",
            )
            println(tex, raw"\midrule")

            base_total = get(get(sdata, ns[1], Dict()), "log_wall_clock_sec", nothing)
            base_abs = get(get(sdata, ns[1], Dict()), "abstraction_time_sec", nothing)
            # If serial data exists, use it as baseline
            if serial_time !== nothing
                base_total = serial_time
            end
            if serial_abs_time !== nothing
                base_abs = serial_abs_time
            end

            # Collect for CSV
            csv_rows = Tuple{Int, Float64, Float64}[]
            for n in ns
                lg = sdata[n]
                total = get(lg, "log_wall_clock_sec", get(lg, "log_timed_sec", NaN))
                abt = get(lg, "abstraction_time_sec", NaN)
                sp_total =
                    (base_total !== nothing && base_total > 0 && total > 0) ?
                    fmt_float(base_total / total) : "—"
                sp_abs =
                    (base_abs !== nothing && base_abs > 0 && abt > 0) ?
                    fmt_float(base_abs / abt) : "—"
                println(
                    tex,
                    "$n & $(isnan(total) ? "—" : fmt_time(total)) & $(isnan(abt) ? "—" : fmt_time(abt)) & $sp_total & $sp_abs \\\\",
                )
                if !isnan(total) && !isnan(abt)
                    push!(csv_rows, (n, total, abt))
                end
            end

            println(tex, raw"\bottomrule")
            println(tex, raw"\end{tabular}")
            println(tex, raw"\end{table}")

            # Write CSV for external plotting
            tag = lowercase(label)
            csv_path = joinpath(report_dir, "scaling_$(tag).csv")
            mkpath(report_dir)
            open(csv_path, "w") do io
                println(io, "n,total_sec,abstraction_sec,speedup_total,speedup_abstraction")
                for (n, total, abt) in csv_rows
                    sp_t =
                        (base_total !== nothing && base_total > 0) ? base_total / total :
                        NaN
                    sp_a = (base_abs !== nothing && base_abs > 0) ? base_abs / abt : NaN
                    println(io, "$n,$total,$abt,$sp_t,$sp_a")
                end
            end

            # Generate gnuplot script and try to produce PNG
            gp_path = joinpath(report_dir, "scaling_$(tag).gp")
            png_path = joinpath(report_dir, "scaling_$(tag).png")
            open(gp_path, "w") do io
                println(
                    io,
                    "set terminal pngcairo size 800,500 enhanced font 'Helvetica,12'",
                )
                println(io, "set output '$(basename(png_path))'")
                println(io, "set title '$label Scaling'")
                println(io, "set xlabel '$xlab'")
                println(io, "set ylabel 'Time (s)'")
                println(io, "set y2label 'Speedup'")
                println(io, "set y2tics")
                println(io, "set key top right")
                println(io, "set grid")
                println(io, "set datafile separator ','")
                println(io, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.2")
                println(io, "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.2")
                println(io, "set style line 3 lc rgb '#00a000' lt 2 lw 2 pt 9 ps 1.2")
                println(
                    io,
                    "plot '$(basename(csv_path))' using 1:2 with linespoints ls 1 title 'Total time', \\",
                )
                println(
                    io,
                    "     '' using 1:3 with linespoints ls 2 title 'Abstraction time', \\",
                )
                return println(
                    io,
                    "     '' using 1:5 axes x1y2 with linespoints ls 3 title 'Speedup (abstr.)'",
                )
            end

            try
                cd(report_dir) do
                    return run(`gnuplot $(basename(gp_path))`)
                end
                if isfile(png_path)
                    println(tex, "\\begin{figure}[H]")
                    println(tex, "\\centering")
                    println(
                        tex,
                        "\\includegraphics[width=0.85\\textwidth]{$(basename(png_path))}",
                    )
                    println(tex, "\\caption{$label scaling: time and speedup vs.\\ $xlab}")
                    println(tex, "\\end{figure}")
                end
            catch
                println("Note: gnuplot not available, skipping $(tag) scaling plot.")
                println("  Data saved to: $csv_path")
            end
        end
    end

    # ── Section 5: Pros & Cons ──
    println(
        tex,
        raw"""
\section{Pros, Cons, and Recommendations}

\subsection{Serial}
\begin{itemize}[leftmargin=2em]
  \item[\textcolor{successgreen}{+}] Simplest code, no synchronization overhead, fully deterministic.
  \item[\textcolor{successgreen}{+}] No worker startup cost.
  \item[\textcolor{red}{--}] Does not exploit multicore hardware at all.
  \item[\textcolor{red}{--}] Abstraction construction time scales linearly with state-space size.
\end{itemize}
\textbf{Best for:} Small state spaces, debugging, baseline measurements.

\subsection{Multithreaded (\texttt{Threads.@threads})}
\begin{itemize}[leftmargin=2em]
  \item[\textcolor{successgreen}{+}] Shared-memory parallelism --- no data serialization cost.
  \item[\textcolor{successgreen}{+}] Easy to use; scales with \texttt{--threads=N}.
  \item[\textcolor{successgreen}{+}] Good speedup on single-node, multicore machines.
  \item[\textcolor{red}{--}] Limited to one machine (shared address space).
  \item[\textcolor{red}{--}] GC stop-the-world pauses affect all threads.
  \item[\textcolor{red}{--}] Thread-safety of user-defined dynamics must be ensured.
\end{itemize}
\textbf{Best for:} Single-node jobs; \texttt{--cpus-per-task} $\geq 4$.

\subsection{Distributed (\texttt{Distributed.jl})}
\begin{itemize}[leftmargin=2em]
  \item[\textcolor{successgreen}{+}] Scales across multiple nodes.
  \item[\textcolor{successgreen}{+}] Each worker has independent GC, memory, fault isolation.
  \item[\textcolor{successgreen}{+}] Natural fit for SLURM \texttt{--ntasks} across nodes.
  \item[\textcolor{red}{--}] Data must be serialized/deserialized for each worker (communication cost).
  \item[\textcolor{red}{--}] Worker startup and \texttt{@everywhere} loading add overhead.
  \item[\textcolor{red}{--}] The final transition merge runs on the master (bottleneck for large abstractions).
\end{itemize}
\textbf{Best for:} Multi-node clusters; large state spaces; embarrassingly parallel abstraction.

\subsection{Hybrid (Distributed + Threaded)}
\begin{itemize}[leftmargin=2em]
  \item[\textcolor{successgreen}{+}] Combines multi-node scaling with intra-node thread parallelism.
  \item[\textcolor{successgreen}{+}] Maximizes hardware utilization ($\text{workers} \times \text{threads}$).
  \item[\textcolor{red}{--}] Most complex configuration.
  \item[\textcolor{red}{--}] Requires careful SLURM tuning: \texttt{--ntasks} $\times$ \texttt{--cpus-per-task}.
  \item[\textcolor{red}{--}] Diminishing returns if communication dominates.
\end{itemize}
\textbf{Best for:} Large-scale HPC; very large state spaces; multi-node with multicore nodes.
""",
    )

    # ── Section 6: SLURM cheat sheet ──
    println(
        tex,
        raw"""
\section{SLURM Configuration Cheat Sheet}

\begin{table}[H]
\centering
\caption{Recommended \texttt{sbatch} parameters per mode}
\begin{tabular}{l l l l l}
\toprule
\textbf{Mode} & \textbf{--nodes} & \textbf{--ntasks} & \textbf{--cpus-per-task} & \textbf{Julia flags} \\
\midrule
Serial      & 1 & 1 & 1          & (none) \\
Threaded    & 1 & 1 & 4--32      & \texttt{--threads=\$SLURM\_CPUS\_PER\_TASK} \\
Distributed & 1--N & 4--32 & 1   & (workers via \texttt{addprocs}) \\
Hybrid      & 1--N & 2--8  & 4--8 & \texttt{--threads=\$SLURM\_CPUS\_PER\_TASK} \\
\bottomrule
\end{tabular}
\end{table}

\textbf{Key rules:}
\begin{itemize}[leftmargin=2em]
  \item \texttt{DIONYSOS\_DISTRIBUTED=true} activates \texttt{Distributed.pmap} partitioned abstraction.
  \item \texttt{DIONYSOS\_THREADED=true} activates \texttt{Threads.@threads} within each process/worker.
  \item \texttt{DIONYSOS\_NPARTS} controls how the state space is partitioned across workers.
  \item Total cores = \texttt{ntasks} $\times$ \texttt{cpus-per-task} $\leq$ node cores $\times$ nodes.
  \item For threaded mode, set \texttt{JULIA\_NUM\_THREADS=\$SLURM\_CPUS\_PER\_TASK}.
  \item Use \texttt{--mem-per-cpu} for distributed/hybrid to scale memory with task count.
\end{itemize}
""",
    )

    # ── Section 7: Scaling guidance ──
    println(
        tex,
        raw"""
\section{Scaling Guidance}

\begin{table}[H]
\centering
\caption{When to use each strategy}
\begin{tabular}{l l l}
\toprule
\textbf{Scenario} & \textbf{Recommended Mode} & \textbf{Rationale} \\
\midrule
Small state space, 1 node      & Serial           & Overhead exceeds benefit \\
Medium state space, 1 node     & Threaded          & Good speedup, minimal complexity \\
Large state space, multi-node  & Distributed       & Distributes memory pressure \\
Very large, multi-node         & Hybrid            & Maximizes all available hardware \\
\bottomrule
\end{tabular}
\end{table}
""",
    )

    # ── End ──
    println(
        tex,
        raw"""
\section{Conclusion}
This benchmark demonstrates the effectiveness of the four parallelism strategies available
in Dionysos for abstraction-based controller synthesis.
The \texttt{UniformGridAbstraction.Optimizer} natively supports serial, threaded, distributed,
and hybrid execution modes via the \texttt{distributed} and \texttt{threaded} parameters.
The optimal strategy depends on the state-space size, available hardware, and the ratio of
computation to communication overhead.

\vfill
\begin{center}
\small
Generated on \today{} by \texttt{generate\_report.jl}\\
Dionysos Parallel Benchmark Suite
\end{center}

\end{document}
""",
    )

    # Write to file
    tex_content = String(take!(tex))
    tex_path = joinpath(report_dir, "benchmark_report.tex")
    mkpath(report_dir)
    open(tex_path, "w") do io
        return write(io, tex_content)
    end
    println("LaTeX report written to: $tex_path")

    # Try to compile
    try
        cd(report_dir) do
            run(`pdflatex -interaction=nonstopmode -halt-on-error benchmark_report.tex`)
            return run(
                `pdflatex -interaction=nonstopmode -halt-on-error benchmark_report.tex`,
            )
        end
        println("PDF report generated: $(joinpath(report_dir, "benchmark_report.pdf"))")
    catch e
        pdf_path = joinpath(report_dir, "benchmark_report.pdf")
        if isfile(pdf_path)
            println("PDF report generated (with warnings): $pdf_path")
        else
            println(
                "pdflatex not available or failed. LaTeX source is ready for manual compilation.",
            )
            println("  cd $report_dir && pdflatex benchmark_report.tex")
        end
    end

    return tex_path
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    script_dir = @__DIR__
    parallel_tests_dir = dirname(script_dir)

    results_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(parallel_tests_dir, "results")
    report_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(parallel_tests_dir, "report")

    println("=" ^ 60)
    println("  DIONYSOS BENCHMARK REPORT GENERATOR")
    println("=" ^ 60)
    println("  Results dir: $results_dir")
    println("  Report dir:  $report_dir")
    println()

    if !isdir(results_dir)
        println("ERROR: Results directory not found: $results_dir")
        println("Run the benchmarks first:")
        println("  bash parallel_tests/run_local.sh <example.jl>")
        exit(1)
    end

    json_files = find_json_files(results_dir)
    if isempty(json_files)
        println("ERROR: No metrics JSON files found in $results_dir")
        println("Make sure the benchmarks were run with the metrics wrapper.")
        exit(1)
    end

    results = load_results(results_dir)
    println("Found results for modes: $(collect(keys(results)))")
    for (mode, runs) in results
        println("  $mode: $(length(runs)) run(s)")
    end
    println()

    log_data = load_log_data(results_dir)
    if !isempty(log_data)
        println("Parsed log files for modes: $(collect(keys(log_data)))")
        for (mode, lg) in log_data
            abt = get(lg, "abstraction_time_sec", nothing)
            abt !== nothing && println("  $mode: abstraction_time = $(abt) s")
        end
        println()
    end

    generate_latex(results, report_dir; log_data = log_data)
    return println("\nDone.")
end

main()
