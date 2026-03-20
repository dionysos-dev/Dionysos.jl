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

    for (regex, key, T) in vcat(patterns, header_patterns, footer_patterns)
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

    # ── Section 4: Results comparison table ──
    println(
        tex,
        raw"""
\section{Performance Comparison}
\subsection{Summary Table}
""",
    )
    println(tex, raw"\begin{table}[H]")
    println(tex, raw"\centering")
    println(tex, raw"\caption{Benchmark Results — all modes}")
    println(tex, raw"\begin{tabular}{l r r r r r r}")
    println(tex, raw"\toprule")
    println(
        tex,
        raw"\textbf{Mode} & \textbf{Total (s)} & \textbf{Speedup} & \textbf{Threads} & \textbf{Workers} & \textbf{Parts} & \textbf{Alloc (MB)} \\\\",
    )
    println(tex, raw"\midrule")

    serial_time = nothing
    if haskey(data, "serial")
        serial_time = get(data["serial"], "elapsed_sec", nothing)
    end

    for m in modes
        haskey(data, m) || continue
        d = data[m]
        t = get(d, "elapsed_sec", 0.0)
        threads = get(d, "julia_threads", 1)
        # Use example_nprocs from captured globals, fall back to log, then to JSON
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

        mode_label = Dict(
            "serial" => "Serial",
            "threaded" => "Threaded",
            "distributed" => "Distributed",
            "hybrid" => "Hybrid",
        )[m]
        println(
            tex,
            "$mode_label & $(fmt_time(t)) & $(speedup)\\texttimes & $threads & $nprocs_actual & $parts & $alloc \\\\",
        )
    end

    println(tex, raw"\bottomrule")
    println(tex, raw"\end{tabular}")
    println(tex, raw"\end{table}")

    # ── Abstraction time comparison (from logs) ──
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

    if !isempty(abs_times)
        println(tex, raw"\subsection{Abstraction Construction Time}")
        println(tex, raw"\begin{table}[H]")
        println(tex, raw"\centering")
        println(
            tex,
            raw"\caption{Abstraction construction time per mode (from example output)}",
        )
        println(tex, raw"\begin{tabular}{l r r r r}")
        println(tex, raw"\toprule")
        println(
            tex,
            raw"\textbf{Mode} & \textbf{Abstraction (s)} & \textbf{Speedup} & \textbf{Threads} & \textbf{Partitions} \\\\",
        )
        println(tex, raw"\midrule")

        for m in modes
            haskey(abs_times, m) || continue
            abt = abs_times[m]
            lg = get(log_data, m, Dict())
            threads = get(lg, "log_threads", 1)
            parts = get(lg, "log_nparts", 1)
            sp = "1.00"
            if serial_abs_time !== nothing && serial_abs_time > 0 && abt > 0
                sp = fmt_float(serial_abs_time / abt)
            end
            mode_label = Dict(
                "serial" => "Serial",
                "threaded" => "Threaded",
                "distributed" => "Distributed",
                "hybrid" => "Hybrid",
            )[m]
            println(
                tex,
                "$mode_label & $(fmt_time(abt)) & $(sp)\\texttimes & $threads & $parts \\\\",
            )
        end

        println(tex, raw"\bottomrule")
        println(tex, raw"\end{tabular}")
        println(tex, raw"\end{table}")
    end

    # ── Section 4: Detailed per-mode analysis ──
    println(tex, raw"\subsection{Detailed Metrics}")

    for m in modes
        haskey(data, m) || continue
        d = data[m]
        mode_label = Dict(
            "serial" => "Serial",
            "threaded" => "Multithreaded",
            "distributed" => "Distributed (Multiprocessing)",
            "hybrid" => "Hybrid (Distributed + Threaded)",
        )[m]

        # Compute speedup for this mode
        speedup_str = "—"
        if serial_time !== nothing && serial_time > 0
            t = get(d, "elapsed_sec", 0.0)
            if t > 0
                speedup_str = fmt_float(serial_time / t)
            end
        end

        println(tex, "\n\\subsubsection{$mode_label}")
        println(tex, raw"\begin{table}[H]")
        println(tex, raw"\centering")
        println(tex, raw"\begin{tabular}{ll}")
        println(tex, raw"\toprule")
        println(tex, raw"\textbf{Metric} & \textbf{Value} \\\\")
        println(tex, raw"\midrule")

        detail_keys = [
            ("elapsed_sec", "Wall-clock time (s)"),
            ("wall_clock_sec", "Total time incl. setup (s)"),
            ("julia_threads", "Julia threads"),
            ("julia_nprocs", "Julia processes"),
            ("julia_nworkers", "Julia workers"),
            ("dionysos_nparts", "Partitions (NPARTS)"),
            ("alloc_MB", "Allocated memory (MB)"),
            ("gc_time_sec", "GC time (s)"),
            ("free_memory_GB", "Free memory (GB)"),
            ("total_memory_GB", "Total memory (GB)"),
            ("hostname", "Hostname"),
            ("pid", "PID"),
        ]

        for (key, label) in detail_keys
            val = get(d, key, nothing)
            val === nothing && continue
            val_str = val isa Bool ? (val ? "Yes" : "No") : escape_latex(string(val))
            println(tex, "$label & \\texttt{$val_str} \\\\")
        end

        println(tex, "Speedup vs.~serial & \\texttt{$(speedup_str)\\texttimes} \\\\")

        # Log-parsed metrics (abstraction time, problem dims, example config)
        lg = get(log_data, m, Dict{String, Any}())
        if !isempty(lg)
            println(tex, raw"\midrule")
            println(tex, raw"\multicolumn{2}{l}{\textbf{From Example Output}} \\\\")
            log_keys = [
                ("abstraction_time_sec", "Abstraction time (s)"),
                ("solve_time_sec", "Abstract problem solve time (s)"),
                ("n_state", "State dimension"),
                ("n_input", "Input dimension"),
                ("n_pos", "Position DoF"),
                ("n_vel", "Velocity DoF"),
            ]
            for (key, label) in log_keys
                val = get(lg, key, nothing)
                val === nothing && continue
                val_str = val isa AbstractFloat ? fmt_time(val) : string(val)
                println(tex, "$label & \\texttt{$(escape_latex(val_str))} \\\\")
            end
            # Abstraction speedup
            if haskey(lg, "abstraction_time_sec") &&
               serial_abs_time !== nothing &&
               serial_abs_time > 0
                abs_sp = fmt_float(serial_abs_time / lg["abstraction_time_sec"])
                println(
                    tex,
                    "Abstraction speedup vs.~serial & \\texttt{$(abs_sp)\\texttimes} \\\\",
                )
            end
        end

        # Captured example globals
        example_keys = [
            ("abstraction_time_sec", "Abstraction time (captured)"),
            ("n_state", "State dim (captured)"),
            ("n_input", "Input dim (captured)"),
            ("example_nprocs", "Example N\\_PROCS"),
            ("example_nparts", "Example N\\_PARTS"),
        ]
        has_example = false
        for (key, label) in example_keys
            val = get(d, key, nothing)
            val === nothing && continue
            if !has_example
                println(tex, raw"\midrule")
                println(
                    tex,
                    raw"\multicolumn{2}{l}{\textbf{Captured Example Globals}} \\\\",
                )
                has_example = true
            end
            val_str = val isa AbstractFloat ? fmt_time(val) : escape_latex(string(val))
            println(tex, "$label & \\texttt{$val_str} \\\\")
        end

        # Worker details
        winfo = get(d, "worker_info", nothing)
        if winfo !== nothing && winfo isa Dict
            println(tex, raw"\midrule")
            println(tex, raw"\multicolumn{2}{l}{\textbf{Worker Details}} \\\\")
            for (wid, info) in winfo
                host = get(info, "host", "?")
                pid = get(info, "pid", "?")
                thr = get(info, "threads", "?")
                println(
                    tex,
                    "Worker $wid & host=\\texttt{$(escape_latex(string(host)))}, pid=$pid, threads=$thr \\\\",
                )
            end
        end

        nodes = get(d, "nodes_used", nothing)
        if nodes !== nothing && nodes isa Vector
            println(tex, "Nodes used & \\texttt{$(escape_latex(join(nodes, ", ")))} \\\\")
        end

        # SLURM info
        has_slurm = false
        for var in [
            "SLURM_JOB_ID",
            "SLURM_NNODES",
            "SLURM_NTASKS",
            "SLURM_CPUS_PER_TASK",
            "SLURM_NODELIST",
        ]
            val = get(d, var, nothing)
            if val !== nothing
                if !has_slurm
                    println(tex, raw"\midrule")
                    println(tex, raw"\multicolumn{2}{l}{\textbf{SLURM Configuration}} \\\\")
                    has_slurm = true
                end
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
