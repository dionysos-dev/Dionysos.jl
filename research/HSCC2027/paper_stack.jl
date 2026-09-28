# Several panels one above the other in a single figure, sharing x and y limits and one x axis.
#
# `latex_axes` in paper_fonts.jl emits a single pgfplots axis, which is right for a figure that
# stands alone. Two figures meant to be *compared* — the primal graph against the dual — read much
# better stacked with the same scales, and a stacked pair also fits a two-column page where two
# separate figures would not.
#
# Include after paper_fonts.jl: this reuses its `crop_png!` and its font preamble.

"Raster `p`'s graphics alone, spanning exactly (x0,x1) x (y0,y1), and return the path written."
function bare_raster(p, path, x0, x1, y0, y1, dpi)
    plot!(
        p;
        framestyle = :box,
        ticks = false,
        xlabel = "",
        ylabel = "",
        title = "",
        xlims = (x0, x1),
        ylims = (y0, y1),
        widen = false,
        left_margin = 0Plots.PlotMeasures.mm,
        right_margin = 0Plots.PlotMeasures.mm,
        top_margin = 0Plots.PlotMeasures.mm,
        bottom_margin = 0Plots.PlotMeasures.mm,
        dpi = dpi,
    )
    savefig(p, path)
    crop_png!(path; pad = 0)          # the image edge is now the frame, i.e. the data limits
    let a = PNGFiles.load(path), k = 2
        PNGFiles.save(path, a[(1 + k):(end - k), (1 + k):(end - k)])
    end
    return path
end

"""
    latex_axes_stacked(ps, out; rows, ...)

`ps` top first. The limits are the union over all panels and every raster is drawn to them, so a
bar at a given height means the same thing in every row. Only the bottom panel carries the x axis
labels; `rows` names each panel inside it, at `row_at` in axis fractions. The legend goes in
`legend_row`, which is worth pointing at whichever panel has room rather than always the first.
"""
function latex_axes_stacked(
    ps,
    out;
    rows = nothing,
    xlab = "x_1",
    ylab = "x_2",
    width_pt = 241,
    panel_pt = 92,
    tick_pt = 9,
    label_pt = 10,
    dpi = 400,
    ylog = false,
    defs = String[],
    legend = nothing,
    legend_pos = "north east",
    legend_row = 1,
    row_at = (0.03, 0.94),
)
    xs, ys = [Plots.xlims(p) for p in ps], [Plots.ylims(p) for p in ps]
    x0, x1 = minimum(first.(xs)), maximum(last.(xs))
    y0, y1 = minimum(first.(ys)), maximum(last.(ys))
    dir = mktempdir()
    for (k, p) in enumerate(ps)
        bare_raster(p, joinpath(dir, "g$k.png"), x0, x1, y0, y1, dpi)
    end

    tickfont = "\\fontsize{$tick_pt}{$(tick_pt + 2)}\\selectfont"
    labelfont = "\\fontsize{$label_pt}{$(label_pt + 2)}\\selectfont"
    lims = "xmin=$x0, xmax=$x1, ymin=$y0, ymax=$y1"

    body = IOBuffer()
    for k in eachindex(ps)
        place = k == 1 ? "" : "at={(ax$(k - 1).below south west)}, anchor=north west,"
        xaxis = k == length(ps) ? "xlabel={\$$xlab\$}," : "xticklabels={},"
        println(body, "\\begin{axis}[")
        println(
            body,
            "  name=ax$k, width=$(width_pt)pt, height=$(panel_pt)pt, scale only axis,",
        )
        isempty(place) || println(body, "  $place")
        println(body, "  axis lines=left, axis line style={line width=0.4pt},")
        println(body, "  enlargelimits=false,")
        ylog && println(body, "  ymode=log,")
        println(body, "  $lims,")
        println(body, "  $xaxis ylabel={\$$ylab\$},")
        println(body, "  tick label style={font=$tickfont},")
        println(body, "  label style={font=$labelfont},")
        println(body, "  tick align=outside,")
        if !isnothing(legend) && k == legend_row
            println(body, "  legend pos=$legend_pos, legend cell align=left,")
            println(
                body,
                "  legend style={font=$tickfont, fill=white, draw=black, rounded corners=0pt},",
            )
        end
        println(body, "]")
        println(body, "\\addplot[forget plot] graphics[$lims]{g$k.png};")
        if !isnothing(legend) && k == legend_row
            for (st, tx) in legend
                println(body, "\\addlegendimage{$st}")
                println(body, "\\addlegendentry{$tx}")
            end
        end
        isnothing(rows) || println(
            body,
            "\\node[anchor=north west, font=$tickfont] at (rel axis cs:$(row_at[1]),$(row_at[2])) {$(rows[k])};",
        )
        println(body, "\\end{axis}")
    end

    open(joinpath(dir, "fig.tex"), "w") do io
        println(io, "\\documentclass[border=1pt]{standalone}")
        println(io, "\\usepackage[T1]{fontenc}")
        println(io, "\\usepackage[tt=false, type1=true]{libertine}")
        println(io, "\\usepackage[libertine]{newtxmath}")
        println(io, "\\usepackage{pgfplots}")
        println(io, "\\pgfplotsset{compat=1.16}")
        println(io, "\\begin{document}")
        foreach(l -> println(io, l), defs)
        println(io, "\\begin{tikzpicture}")
        print(io, String(take!(body)))
        println(io, "\\end{tikzpicture}")
        return println(io, "\\end{document}")
    end

    cd(dir) do
        Base.run(`pdflatex -interaction=nonstopmode -halt-on-error fig.tex`)
        return nothing
    end
    Base.run(
        `pdftocairo -png -r $dpi -singlefile $(joinpath(dir, "fig.pdf")) $(splitext(out)[1])`,
    )
    crop_png!(out)
    rm(dir; recursive = true, force = true)
    return out
end
