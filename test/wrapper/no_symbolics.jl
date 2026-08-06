module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

# The point of moving the front-end out of the extension is that it works with no optional
# dependency loaded. That cannot be checked from inside this session: `runtests.jl` includes
# every file into a single process, and several of them do `using Symbolics`, which activates
# the extension for everyone afterwards. So the check runs in a fresh subprocess.
#
# The subprocess must inherit the *active* environment, not `test/`: `Pkg.test()` runs in a
# temporary sandbox, where `test/Project.toml` is not instantiated and `using Dionysos` fails.
const PROJECT =
    something(Base.active_project(), joinpath(dirname(dirname(pathof(Dionysos))), "test"))

# Deliberately loads neither Symbolics nor MathOptSymbolicAD.
const SCRIPT = """
using Dionysos, JuMP, StaticArrays
import MathOptInterface as MOI
const MP = Dionysos.Mapping
const AB = Dionysos.Optim.Abstraction

if Base.get_extension(Dionysos, :DionysosMathOptSymbolicADExt) !== nothing
    error("the MathOptSymbolicAD extension is loaded; this check proves nothing")
end

model = direct_model(Dionysos.Optimizer())
@variable(model, -1.0 <= x <= 1.0, start = -0.75)
@variable(model, -1.0 <= u <= 1.0)
@constraint(model, d(x) == u)
@constraint(model, fin(x) in MOI.Interval(-0.5, 0.5))

set_attribute(model, "dynamics_backend", Dionysos.Wrapper.NonlinearEvaluatorBackend())
set_attribute(model, "time_step", 1.0)
set_attribute(model, "approx_mode", AB.UniformGridAbstraction.CENTER_SIMULATION)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.25)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.5)))
set_attribute(model, "print_level", 0)

optimize!(model)
is_solved_and_feasible(model) || error("no controller")
Dionysos.simulate(model, SVector(-0.75); nsteps = 20)
println("WRAPPER-OK")
"""

@testset "the front-end works with no optional dependency loaded" begin
    # `∂` and `final` are written as `d`/`fin` in the script above only to keep the command
    # line free of non-ASCII quoting surprises across platforms.
    script = replace(SCRIPT, "d(x)" => "∂(x)", "fin(x)" => "final(x)")
    command = `$(Base.julia_cmd()) --project=$(PROJECT) -e $script`

    buffer = IOBuffer()
    ok = success(pipeline(command; stdout = buffer, stderr = buffer))
    output = String(take!(buffer))
    ok || println(output)

    @test ok
    @test occursin("WRAPPER-OK", output)
end

end # module TestMain
