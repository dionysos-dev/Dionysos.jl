# Generic asset runner: includes an example driver and saves its animation and
# figure without modifying the driver. Usage:
#   julia --project=test presentation/run_driver.jl <driver.jl> <out.gif>
#
# The animation is found, in order: the include's return value (most drivers end
# with the `animate_trajectory_dashboard` call), an already-rendered AnimatedGif
# (its temp file is copied), or any `Plots.Animation` bound in Main.

ENV["GKSwstype"] = "100"

driver = abspath(ARGS[1])
out = abspath(ARGS[2])

result = include(driver)

import Plots

function save_animation(x, out)
    if x isa Plots.Animation
        Plots.gif(x, out; fps = 8)
        return true
    elseif x isa Plots.AnimatedGif
        cp(x.filename, out; force = true)
        return true
    end
    return false
end

saved = save_animation(result, out)
if !saved
    for name in names(Main; all = true)
        isdefined(Main, name) || continue
        (name in (:result, :saved)) && continue
        if save_animation(getfield(Main, name), out)
            saved = true
            break
        end
    end
end
println(saved ? "saved $out" : "no animation found")

if @isdefined(fig)
    png_out = replace(out, r"\.gif$" => ".png")
    Plots.savefig(fig, png_out)
    println("saved ", png_out)
end
