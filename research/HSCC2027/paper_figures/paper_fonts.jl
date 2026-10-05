# Put the paper's own type on a figure whose graphics are too heavy for LaTeX to draw.
#
# PGFPlotsX gives exact acmart fonts by emitting the whole figure as LaTeX source, which works for a
# dozen polygons and not for 22 438 — 20 000 took lualatex past 45 minutes with no page. The way to
# have both is to let GR raster the graphics with no axes at all, then have LaTeX draw only the
# frame, the ticks and the labels around that image. LaTeX sees one `\includegraphics`, so it
# compiles in seconds, and every glyph is Libertine.
#
# Needs `pdflatex` and `pdftocairo` on PATH. `has_latex()` says whether they are there, so a caller
# can fall back to the plain GR figure.

import PNGFiles

has_latex() = !isnothing(Sys.which("pdflatex")) && !isnothing(Sys.which("pdftocairo"))

"Trim the uniform border off a PNG. `pad` leaves a hair of margin so nothing is shaved."
function crop_png!(path; pad = 3)
    img = PNGFiles.load(path)
    bg = img[1, 1]
    rows = findall(i -> any(!=(bg), view(img, i, :)), axes(img, 1))
    cols = findall(j -> any(!=(bg), view(img, :, j)), axes(img, 2))
    (isempty(rows) || isempty(cols)) && return nothing
    r = max(first(rows) - pad, 1):min(last(rows) + pad, size(img, 1))
    c = max(first(cols) - pad, 1):min(last(cols) + pad, size(img, 2))
    PNGFiles.save(path, img[r, c])
    return nothing
end

"""
    latex_axes(p, out; xlab, ylab, width_pt, tick_pt, label_pt, dpi)

Draw `p`'s graphics as a bare raster, then wrap it in a pgfplots axis typeset by `pdflatex` in
acmart's fonts, and leave the result at `out` as a single cropped PNG.

The axis is given the figure's own limits and `enlargelimits=false`, so the image lines up with the
ticks exactly; `p` is stripped of its frame first, or GR's axes would be baked into the raster.
"""
function latex_axes(
    p,
    out;
    xlab = "x_1",
    ylab = "x_2",
    width_pt = 241,
    height_pt = nothing,
    tick_pt = 11,
    label_pt = 14,
    dpi = 400,
    legend = nothing,
    legend_pos = "north east",   # the disc leaves the corners empty
    ylog = false,
    # Plots colour names mean nothing to TikZ; pass `\definecolor` lines for any the legend uses
    # beyond xcolor's own (black, blue, red, ...).
    defs = String[],
)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    h = something(height_pt, round(Int, width_pt * (y1 - y0) / (x1 - x0)))
    dir = mktempdir()
    img = joinpath(dir, "graphics.png")

    # The raster must span the data limits EXACTLY, or the cells sit off the ticks. Chasing GR's
    # padding with negative margins does not do it. Instead the figure is drawn with a box frame at
    # the limits and nothing else, the PNG is cropped to that frame — which is then by construction
    # the outermost ink — and the frame line itself is shaved off.
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
    savefig(p, img)
    crop_png!(img; pad = 0)          # now the image edge is the frame, i.e. the data limits
    let a = PNGFiles.load(img), k = 2
        PNGFiles.save(img, a[(1 + k):(end - k), (1 + k):(end - k)])
    end

    tex = joinpath(dir, "fig.tex")
    open(tex, "w") do io
        return print(
            io,
            """
\\documentclass[border=1pt]{standalone}
\\usepackage[T1]{fontenc}
\\usepackage[tt=false, type1=true]{libertine}
\\usepackage[libertine]{newtxmath}
\\usepackage{pgfplots}
\\pgfplotsset{compat=1.16}
\\begin{document}
$(isempty(defs) ? "" : join(defs, "\n") * "\n")\\begin{tikzpicture}
\\begin{axis}[
  width=$(width_pt)pt, height=$(h)pt, scale only axis,
  axis lines=left, axis line style={line width=0.4pt},
  enlargelimits=false,
$(ylog ? "  ymode=log,\n" : "")  xmin=$x0, xmax=$x1, ymin=$y0, ymax=$y1,
  xlabel={\$$xlab\$}, ylabel={\$$ylab\$},
  tick label style={font=\\fontsize{$tick_pt}{$(tick_pt + 2)}\\selectfont},
  label style={font=\\fontsize{$label_pt}{$(label_pt + 2)}\\selectfont},
  tick align=outside,
$(isnothing(legend) ? "" : "  legend pos=$legend_pos, legend cell align=left,
  legend style={font=\\fontsize{$tick_pt}{$(tick_pt + 2)}\\selectfont, fill=white, draw=black, rounded corners=0pt},
")]
\\addplot[forget plot] graphics[xmin=$x0, xmax=$x1, ymin=$y0, ymax=$y1]{graphics.png};
$(isnothing(legend) ? "" : join(["\\addlegendimage{$st}
\\addlegendentry{$tx}" for (st, tx) in legend], "
") * "
")
\\end{axis}
\\end{tikzpicture}
\\end{document}
""",
        )
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
