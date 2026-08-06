# ----------------------------
# Plotting recipes
# ----------------------------

@recipe function f(
    problem::AlternatingSimulationProblem;
    domain_color = :gray,
    region_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end
    @series begin
        label := "Region"
        color := region_color
        problem.state_set
    end
end

@recipe function f(
    problem::BisimulationQuotientProblem;
    region_color = :gray,
    region_alpha = 0.15,
    observation_region_alpha = 0.4,
    observation_colors = [:red, :green, :orange, :purple, :brown],
    plot_region = true,
)
    if plot_region
        @series begin
            label := "Region"
            color := region_color
            fillalpha := region_alpha
            problem.state_set
        end
    end
    for (i, R) in enumerate(problem.observation_regions)
        @series begin
            label := "O $i"
            color := observation_colors[mod1(i, length(observation_colors))]
            fillalpha := observation_region_alpha
            R
        end
    end
end

@recipe function f(
    problem::OptimalControlProblem;
    domain_color = :gray,
    initial_set_color = :green,
    target_set_color = :red,
    safe_set_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end
    if problem.safe_set !== nothing
        @series begin
            label := "Safe set"
            color := safe_set_color
            problem.safe_set
        end
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
    @series begin
        label := "Target set"
        color := target_set_color
        problem.target_set
    end
end

@recipe function f(
    problem::SafetyProblem;
    domain_color = :gray,
    initial_set_color = :green,
    safe_set_color = :lightgray,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end
    @series begin
        label := "Safe set"
        color := safe_set_color
        problem.safe_set
    end
    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
end

@recipe function f(
    problem::ReachAndStayProblem;
    domain_color = :gray,
    safe_set_color = :lightgray,
    target_set_color = :red,
    initial_set_color = :green,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end

    @series begin
        label := "Safe set"
        color := safe_set_color
        problem.safe_set
    end

    @series begin
        label := "Target set"
        color := target_set_color
        problem.target_set
    end

    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end
end

@recipe function f(
    problem::CoSafeLTLProblem;
    domain_color = :gray,
    initial_set_color = :green,
    ap_colors = Dict{Symbol, Any}(),
    obs_color = :red,
)
    @series begin
        label := "Domain"
        color := domain_color
        MS.stateset(problem.system)
    end

    @series begin
        label := "Initial set"
        color := initial_set_color
        problem.initial_set
    end

    for (ap, setX) in problem.labeling
        color_ap = haskey(ap_colors, ap) ? ap_colors[ap] : :blue
        @series begin
            label := String(ap)
            color := color_ap
            setX
        end
    end
end
