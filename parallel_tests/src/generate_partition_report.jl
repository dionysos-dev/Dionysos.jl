#!/usr/bin/env julia
"""
Generate a LaTeX/PDF report for the partitioned abstraction pipeline.

Usage:
    julia --project=<PROJECT> generate_partition_report.jl <partitions_dir> [report_dir]

Reads `partition_*.jld2` files from `partitions_dir` and, when available,
`assembly_metadata.json`. Produces a LaTeX report summarizing per-partition
performance, estimated total time, and parallel gain.
"""

using Dates
using JLD2
using JSON
using Printf
using Statistics

safe_get(d::Dict, key::String, default = nothing) = get(d, key, default)

function fmt_float(x; digits = 2)
    x isa Number || return string(x)
    return @sprintf("%.*f", digits, Float64(x))
end

function fmt_int(x)
    x isa Number || return string(x)
    return string(round(Int, x))
end

function fmt_time(x)
    x isa Number || return string(x)
    val = Float64(x)
    if val >= 3600
        h = floor(Int, val / 3600)
        m = floor(Int, (val - 3600h) / 60)
        s = val - 3600h - 60m
        return @sprintf("%dh %02dm %05.2fs", h, m, s)
    elseif val >= 60
        m = floor(Int, val / 60)
        s = val - 60m
        return @sprintf("%dm %05.2fs", m, s)
    else
        return @sprintf("%.2f s", val)
    end
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

function partition_index(path::String)
    name = basename(path)
    m = match(r"^partition_(\d+)\.jld2$", name)
    m === nothing && return nothing
    return parse(Int, m[1])
end

function find_partition_files(partitions_dir::String)
    files = String[]
    for entry in readdir(partitions_dir; join = true)
        idx = partition_index(entry)
        idx === nothing && continue
        push!(files, entry)
    end
    sort!(files; by = x -> partition_index(x))
    return files
end

function load_partition_data(partitions_dir::String)
    rows = Dict{String, Any}[]
    for path in find_partition_files(partitions_dir)
        data = load(path)
        meta = data["metadata"]::Dict{String, Any}
        push!(
            rows,
            Dict{String, Any}(
                "path" => path,
                "partition_idx" => safe_get(meta, "partition_idx", partition_index(path)),
                "n_source_states" => safe_get(meta, "n_source_states", 0),
                "n_transitions" => safe_get(meta, "n_transitions", 0),
                "elapsed_compute_sec" => safe_get(meta, "elapsed_compute_sec", 0.0),
                "setup_time_sec" => safe_get(meta, "setup_time_sec", 0.0),
                "partition_time_sec" => safe_get(meta, "partition_time_sec", 0.0),
                "rbd_time_sec" => safe_get(meta, "rbd_time_sec", 0.0),
                "rbd_call_count" => safe_get(meta, "rbd_call_count", 0),
                "total_wall_clock_sec" => safe_get(
                    meta,
                    "total_wall_clock_sec",
                    safe_get(meta, "elapsed_compute_sec", 0.0),
                ),
                "hostname" => safe_get(meta, "hostname", "unknown"),
                "julia_threads" => safe_get(meta, "julia_threads", "unknown"),
                "timestamp" => safe_get(meta, "timestamp", "unknown"),
                "strategy" => safe_get(meta, "strategy", "unknown"),
                "setup_script" => safe_get(meta, "setup_script", "unknown"),
            ),
        )
    end
    return rows
end

function load_assembly_metadata(partitions_dir::String)
    path = joinpath(partitions_dir, "assembly_metadata.json")
    return isfile(path) ? JSON.parsefile(path) : nothing
end

function mean_or_zero(values)
    isempty(values) && return 0.0
    return mean(Float64.(values))
end

function max_or_zero(values)
    isempty(values) && return 0.0
    return maximum(Float64.(values))
end

function min_or_zero(values)
    isempty(values) && return 0.0
    return minimum(Float64.(values))
end

function sum_or_zero(values)
    isempty(values) && return 0.0
    return sum(Float64.(values))
end

function summarize(rows::Vector{Dict{String, Any}}, assembly_meta)
    nparts_expected = if assembly_meta !== nothing
        Int(safe_get(assembly_meta, "nparts", length(rows)))
    elseif isempty(rows)
        0
    else
        maximum(Int(row["partition_idx"]) for row in rows)
    end

    completed_indices = sort([Int(row["partition_idx"]) for row in rows])
    expected_indices = collect(1:nparts_expected)
    missing_indices =
        assembly_meta !== nothing ? Int.(safe_get(assembly_meta, "missing_parts", Int[])) :
        setdiff(expected_indices, completed_indices)

    compute_secs = [Float64(row["elapsed_compute_sec"]) for row in rows]
    wall_secs = [Float64(row["total_wall_clock_sec"]) for row in rows]
    setup_secs = [Float64(row["setup_time_sec"]) for row in rows]
    partition_secs = [Float64(row["partition_time_sec"]) for row in rows]
    rbd_secs = [Float64(row["rbd_time_sec"]) for row in rows]
    source_states = [Float64(row["n_source_states"]) for row in rows]
    transitions = [Float64(row["n_transitions"]) for row in rows]

    sum_compute = sum_or_zero(compute_secs)
    max_compute = max_or_zero(compute_secs)
    max_wall = max_or_zero(wall_secs)
    assembly_total =
        assembly_meta === nothing ? 0.0 :
        Float64(safe_get(assembly_meta, "total_wall_clock_sec", 0.0))
    assembly_add =
        assembly_meta === nothing ? 0.0 :
        Float64(safe_get(assembly_meta, "add_transitions_sec", 0.0))
    estimated_partition_makespan = max_wall > 0 ? max_wall : max_compute
    estimated_total_wall = estimated_partition_makespan + assembly_total
    partition_gain =
        estimated_partition_makespan > 0 ? sum_compute / estimated_partition_makespan : 0.0
    total_gain = estimated_total_wall > 0 ? sum_compute / estimated_total_wall : 0.0

    return Dict{String, Any}(
        "nparts_expected" => nparts_expected,
        "nparts_completed" => length(rows),
        "missing_indices" => missing_indices,
        "sum_compute_sec" => sum_compute,
        "avg_compute_sec" => mean_or_zero(compute_secs),
        "min_compute_sec" => min_or_zero(compute_secs),
        "max_compute_sec" => max_compute,
        "avg_wall_sec" => mean_or_zero(wall_secs),
        "max_wall_sec" => max_wall,
        "avg_setup_sec" => mean_or_zero(setup_secs),
        "avg_partition_sec" => mean_or_zero(partition_secs),
        "avg_rbd_sec" => mean_or_zero(rbd_secs),
        "sum_rbd_sec" => sum_or_zero(rbd_secs),
        "avg_source_states" => mean_or_zero(source_states),
        "avg_transitions" => mean_or_zero(transitions),
        "sum_source_states" => sum_or_zero(source_states),
        "sum_transitions" => sum_or_zero(transitions),
        "estimated_partition_makespan_sec" => estimated_partition_makespan,
        "assembly_total_wall_sec" => assembly_total,
        "assembly_add_sec" => assembly_add,
        "estimated_total_wall_sec" => estimated_total_wall,
        "partition_parallel_gain" => partition_gain,
        "end_to_end_parallel_gain" => total_gain,
        "parallel_efficiency" => length(rows) > 0 ? partition_gain / length(rows) : 0.0,
    )
end

function generate_latex(
    rows::Vector{Dict{String, Any}},
    summary::Dict{String, Any},
    report_dir::String;
    assembly_meta = nothing,
)
    isempty(rows) && error("No partition results found.")

    sample = rows[1]
    setup_script = escape_latex(string(safe_get(sample, "setup_script", "unknown")))
    strategy = escape_latex(string(safe_get(sample, "strategy", "unknown")))
    report_title = escape_latex(basename(dirname(String(sample["path"]))))
    assembly_present = assembly_meta !== nothing

    tex = IOBuffer()
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
\usepackage{array}
\usepackage{float}
\usepackage{fancyhdr}

\setlength{\headheight}{14pt}
\addtolength{\topmargin}{-2pt}

\definecolor{headerblue}{RGB}{41,65,122}
\definecolor{lightgray}{RGB}{240,240,240}

\pagestyle{fancy}
\fancyhf{}
\rhead{Dionysos Partition Pipeline}
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
        "\\title{\\textcolor{headerblue}{\\textbf{Partition Pipeline Report}}\\\\[0.5em]",
    )
    println(tex, "       \\large \\texttt{$report_title}}")
    println(tex, "\\author{Auto-generated by \\texttt{generate\\_partition\\_report.jl}}")
    println(tex, "\\date{", escape_latex(string(Dates.now())), "}")
    println(tex, raw"\begin{document}")
    println(tex, raw"\maketitle")
    println(tex, raw"\tableofcontents")
    println(tex, raw"\newpage")

    println(tex, raw"\section{Executive Summary}")
    println(
        tex,
        "This report summarizes a partitioned Dionysos abstraction run with ",
        "$(fmt_int(summary["nparts_completed"])) completed partitions out of ",
        "$(fmt_int(summary["nparts_expected"])). ",
        "Partition-level compute time sums to $(fmt_time(summary["sum_compute_sec"])), ",
        "while the estimated partition makespan is $(fmt_time(summary["estimated_partition_makespan_sec"])). ",
        "The resulting partition parallel gain estimate is $(fmt_float(summary["partition_parallel_gain"]))\\texttimes.",
    )
    if assembly_present
        println(
            tex,
            " Including assembly, the estimated end-to-end wall time is ",
            "$(fmt_time(summary["estimated_total_wall_sec"])) and the end-to-end gain estimate is ",
            "$(fmt_float(summary["end_to_end_parallel_gain"]))\\texttimes.",
        )
    end
    println(
        tex,
        raw"""

\paragraph{Interpretation.}
The \emph{partition parallel gain} is estimated as
\[
\frac{\sum_i t_i^{(\mathrm{compute})}}{\max_i t_i^{(\mathrm{wall})}}
\]
using per-partition wall-clock when available and per-partition compute time otherwise.
This is an estimate of useful parallelism in the partition stage, not a scheduler-aware SLURM makespan.
""",
    )

    println(tex, raw"\section{Run Configuration}")
    println(tex, raw"\begin{table}[H]")
    println(tex, raw"\centering")
    println(tex, raw"\begin{tabular}{ll}")
    println(tex, raw"\toprule")
    println(tex, "\\textbf{Property} & \\textbf{Value} \\\\")
    println(tex, raw"\midrule")
    println(tex, "Setup script & \\texttt{$setup_script} \\\\")
    println(tex, "Strategy & \\texttt{$strategy} \\\\")
    println(
        tex,
        "Completed partitions & $(fmt_int(summary["nparts_completed"])) / $(fmt_int(summary["nparts_expected"])) \\\\",
    )
    println(
        tex,
        "Julia threads per task & $(escape_latex(string(safe_get(sample, "julia_threads", "unknown")))) \\\\",
    )
    println(
        tex,
        "Representative host & \\texttt{$(escape_latex(string(safe_get(sample, "hostname", "unknown"))))} \\\\",
    )
    if assembly_present
        println(
            tex,
            "Assembly metadata & \\texttt{$(escape_latex(string(safe_get(assembly_meta, "output_dir", "present"))))} \\\\",
        )
    else
        println(tex, "Assembly metadata & not found \\\\")
    end
    println(tex, raw"\bottomrule")
    println(tex, raw"\end{tabular}")
    println(tex, raw"\end{table}")

    println(tex, raw"\section{Aggregate Metrics}")
    println(tex, raw"\begin{table}[H]")
    println(tex, raw"\centering")
    println(tex, raw"\begin{tabular}{lll}")
    println(tex, raw"\toprule")
    println(tex, "\\textbf{Metric} & \\textbf{Value} & \\textbf{Comment} \\\\")
    println(tex, raw"\midrule")
    println(
        tex,
        "Average partition compute time & $(fmt_time(summary["avg_compute_sec"])) & Mean of \\texttt{elapsed\\_compute\\_sec} across completed partitions \\\\",
    )
    println(
        tex,
        "Average partition wall time & $(fmt_time(summary["avg_wall_sec"])) & Mean of partition wall-clock metadata \\\\",
    )
    println(
        tex,
        "Average setup time & $(fmt_time(summary["avg_setup_sec"])) & Per-partition setup cost before transition collection \\\\",
    )
    println(
        tex,
        "Average partitioning time & $(fmt_time(summary["avg_partition_sec"])) & Time spent slicing the source-state space \\\\",
    )
    println(
        tex,
        "Average RBD simulation time & $(fmt_time(summary["avg_rbd_sec"])) & Mean rigid-body simulation time per partition \\\\",
    )
    println(
        tex,
        "Average source states per partition & $(fmt_float(summary["avg_source_states"]; digits = 1)) & Mean local source-state count \\\\",
    )
    println(
        tex,
        "Average transitions per partition & $(fmt_float(summary["avg_transitions"]; digits = 1)) & Mean local transition count \\\\",
    )
    println(
        tex,
        "Total partition compute time & $(fmt_time(summary["sum_compute_sec"])) & Serial-equivalent sum of partition compute times \\\\",
    )
    println(
        tex,
        "Estimated partition makespan & $(fmt_time(summary["estimated_partition_makespan_sec"])) & Maximum partition wall-clock or compute time \\\\",
    )
    println(
        tex,
        "Partition parallel gain & $(fmt_float(summary["partition_parallel_gain"]))\\texttimes & Sum of compute times divided by estimated makespan \\\\",
    )
    println(
        tex,
        "Parallel efficiency & $(fmt_float(100 * summary["parallel_efficiency"]))\\% & Gain divided by completed partition count \\\\",
    )
    if assembly_present
        println(
            tex,
            "Assembly total wall time & $(fmt_time(summary["assembly_total_wall_sec"])) & From \\texttt{assembly\\_metadata.json} \\\\",
        )
        println(
            tex,
            "Assembly add-transitions time & $(fmt_time(summary["assembly_add_sec"])) & Time spent in \\texttt{add\\_transitions!} \\\\",
        )
        println(
            tex,
            "Estimated end-to-end wall time & $(fmt_time(summary["estimated_total_wall_sec"])) & Estimated partition makespan plus assembly wall time \\\\",
        )
        println(
            tex,
            "End-to-end parallel gain & $(fmt_float(summary["end_to_end_parallel_gain"]))\\texttimes & Sum of partition compute times divided by estimated end-to-end time \\\\",
        )
    end
    println(tex, raw"\bottomrule")
    println(tex, raw"\end{tabular}")
    println(tex, raw"\end{table}")

    if !isempty(summary["missing_indices"])
        println(tex, raw"\section{Missing Partitions}")
        println(
            tex,
            "The following partition indices were expected but not found: \\texttt{",
            escape_latex(join(string.(summary["missing_indices"]), ", ")),
            "}.",
        )
    end

    println(tex, raw"\section{Per-Partition Results}")
    println(tex, raw"\rowcolors{2}{lightgray}{white}")
    println(tex, raw"\begin{longtable}{rrrrrrr}")
    println(tex, "\\caption{Per-partition metrics}\\\\")
    println(tex, raw"\toprule")
    println(
        tex,
        "\\textbf{Part} & \\textbf{States} & \\textbf{Transitions} & \\textbf{Compute (s)} & \\textbf{Wall (s)} & \\textbf{RBD (s)} & \\textbf{Calls} \\\\",
    )
    println(tex, raw"\midrule")
    println(tex, raw"\endfirsthead")
    println(tex, raw"\toprule")
    println(
        tex,
        "\\textbf{Part} & \\textbf{States} & \\textbf{Transitions} & \\textbf{Compute (s)} & \\textbf{Wall (s)} & \\textbf{RBD (s)} & \\textbf{Calls} \\\\",
    )
    println(tex, raw"\midrule")
    println(tex, raw"\endhead")
    println(tex, raw"\bottomrule")
    println(tex, raw"\endfoot")
    for row in sort(rows; by = r -> Int(r["partition_idx"]))
        println(
            tex,
            "$(fmt_int(row["partition_idx"])) & $(fmt_int(row["n_source_states"])) & " *
            "$(fmt_int(row["n_transitions"])) & $(fmt_float(row["elapsed_compute_sec"])) & " *
            "$(fmt_float(row["total_wall_clock_sec"])) & $(fmt_float(row["rbd_time_sec"])) & " *
            "$(fmt_int(row["rbd_call_count"])) \\\\",
        )
    end
    println(tex, raw"\end{longtable}")

    println(tex, raw"\end{document}")

    tex_path = joinpath(report_dir, "partition_report.tex")
    mkpath(report_dir)
    open(tex_path, "w") do io
        return write(io, String(take!(tex)))
    end
    println("LaTeX report written to: $tex_path")

    try
        cd(report_dir) do
            run(`pdflatex -interaction=nonstopmode -halt-on-error partition_report.tex`)
            return run(
                `pdflatex -interaction=nonstopmode -halt-on-error partition_report.tex`,
            )
        end
        println("PDF report generated: $(joinpath(report_dir, "partition_report.pdf"))")
    catch
        pdf_path = joinpath(report_dir, "partition_report.pdf")
        if isfile(pdf_path)
            println("PDF report generated (with warnings): $pdf_path")
        else
            println(
                "pdflatex not available or failed. LaTeX source is ready for manual compilation.",
            )
            println("  cd $report_dir && pdflatex partition_report.tex")
        end
    end

    return tex_path
end

function main()
    parallel_tests_dir = dirname(@__DIR__)
    partitions_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(parallel_tests_dir, "results")
    report_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(partitions_dir, "report")

    println("=" ^ 60)
    println("  DIONYSOS PARTITION REPORT GENERATOR")
    println("=" ^ 60)
    println("  Partitions dir: $partitions_dir")
    println("  Report dir:     $report_dir")
    println()

    isdir(partitions_dir) || error("Partitions directory not found: $partitions_dir")

    rows = load_partition_data(partitions_dir)
    isempty(rows) && error("No partition_*.jld2 files found in $partitions_dir")

    assembly_meta = load_assembly_metadata(partitions_dir)
    summary = summarize(rows, assembly_meta)

    println(
        "Completed partitions: $(summary["nparts_completed"]) / $(summary["nparts_expected"])",
    )
    println(
        "Estimated partition makespan: $(fmt_time(summary["estimated_partition_makespan_sec"]))",
    )
    println("Partition parallel gain: $(fmt_float(summary["partition_parallel_gain"]))x")
    assembly_meta !== nothing && println("Assembly metadata found.")
    println()

    generate_latex(rows, summary, report_dir; assembly_meta = assembly_meta)
    return println("\nDone.")
end

main()
