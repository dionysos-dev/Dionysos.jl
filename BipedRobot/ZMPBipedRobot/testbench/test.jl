
stepR = sf.stepR
stepL = sf.stepL

stepR_plot = reduce(hcat, stepR)
stepL_plot = reduce(hcat, stepL)
step_plot = reduce(vcat, [stepL_plot, stepR_plot])

plt = plot( title = "Swing Foot Trajectory",
            xlabel = "X[m]", ylabel = "Z[m]",
            #layout = (1,2), 
            dpi = 600)
i = 1 
for foot in ["left", "right"]
    plot!(step_plot[i + 2 * (i - 1), :], step_plot[i + 1 + 2 * (i - 1), :], step_plot[i + 2 + 2 * (i - 1), :], label = foot)
    i = i + 1; 
end 
plt