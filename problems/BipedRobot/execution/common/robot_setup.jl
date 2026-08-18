const BIPED_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const PROJECT_DIR = BIPED_ROOT
const DEPS_DIR = joinpath(BIPED_ROOT, "deps")
const OUT_DIR = joinpath(BIPED_ROOT, "out")

function selected_robot_model()
    return get(ENV, "DIONYSOS_ROBOT_MODEL", "6D")
end

function selected_robot_problem_path()
    model = selected_robot_model()

    if model == "4D"
        return joinpath(BIPED_ROOT, "4D_model", "robot_problem.jl")
    elseif model == "6D"
        return joinpath(BIPED_ROOT, "6D_model", "robot_problem.jl")
    elseif model == "8D"
        return joinpath(BIPED_ROOT, "8D_model", "robot_problem.jl")
    else
        error("Unknown DIONYSOS_ROBOT_MODEL=$model. Expected 4D, 6D or 8D.")
    end
end

function selected_robot_urdf()
    return joinpath(DEPS_DIR, "ZMP_2DBipedRobot_nodamping.urdf")
end

function default_transition_outdir()
    model = selected_robot_model()
    return joinpath(OUT_DIR, model, "transitions")
end

function default_abstraction_file()
    model = selected_robot_model()
    return joinpath(OUT_DIR, model, "robot_abstraction.jld2")
end
