using Plots

# ------------------------------------------------
# Data
# ------------------------------------------------
nprocs_vals = [1, 2, 3, 10, 20, 40, 60]

worker_creation = [6.1, 5.8, 10.0, 8.5, 8.9, 11.0, 8.9]
worker_package_load = [6.0, 6.0, 6.0, 6.0, 8.3, 12.0, 8.8]
warmup = [24.4, 24.4, 25.0, 31.0, 37.0, 58.0, 65.9]
abstraction_construction = [9121.0, 4535.0, 2900.0, 949.0, 510.0, 276.0, 180.0]
total = [9149.0, 4563.0, 2930.0, 984.0, 551.0, 338.0, 250.0]

ideal = total[1] .* (nprocs_vals[1] ./ nprocs_vals)

# ------------------------------------------------
# Figure 1 — linear scale
# ------------------------------------------------
p1 = plot(;
    xlabel = "Number of processes",
    ylabel = "Time (s)",
    title = "Distributed abstraction timing breakdown",
    legend = :topright,
    grid = true,
    lw = 2,
)

plot!(p1, nprocs_vals, worker_creation; label = "worker creation", marker = :circle)
plot!(p1, nprocs_vals, worker_package_load; label = "worker package load", marker = :square)
plot!(p1, nprocs_vals, warmup; label = "warmup", marker = :diamond)
plot!(
    p1,
    nprocs_vals,
    abstraction_construction;
    label = "abstraction construction",
    marker = :utriangle,
)
plot!(p1, nprocs_vals, total; label = "total", marker = :star5)
plot!(p1, nprocs_vals, ideal; label = "ideal 1/N", linestyle = :dash, marker = :none)

display(p1)
# savefig(p1, "timing_linear.png")

# ------------------------------------------------
# Figure 2 — log-log scale
# ------------------------------------------------
p2 = plot(;
    xlabel = "Number of processes",
    ylabel = "Time (s)",
    title = "Scaling behavior (log₂-log₂)",
    legend = :topright,
    grid = true,
    lw = 2,
    xscale = :log2,
    yscale = :log2,
)

plot!(p2, nprocs_vals, worker_creation; label = "worker creation", marker = :circle)
plot!(p2, nprocs_vals, worker_package_load; label = "worker package load", marker = :square)
plot!(p2, nprocs_vals, warmup; label = "warmup", marker = :diamond)
plot!(
    p2,
    nprocs_vals,
    abstraction_construction;
    label = "abstraction construction",
    marker = :utriangle,
)
plot!(p2, nprocs_vals, total; label = "total", marker = :star5)
plot!(p2, nprocs_vals, ideal; label = "ideal 1/N", linestyle = :dash, marker = :none)
# savefig(p2, "timing_loglog.png")
display(p2)
