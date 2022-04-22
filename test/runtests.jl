using Dionysos

include("./utils/data_structures/BDD/test_BDD.jl")
include("./utils/data_structures/BDD/test_inttupleset.jl")
# See https://github.com/jump-dev/SDPA.jl/issues/34
# Once this issue is resolved, this if must be removed
if VERSION == v"1.6.3"
    include("./examples/test_state_feedback_trans.jl")
end
include("./examples/test_gol_lazar_belta.jl")
include("./domain/test_griddomain.jl")
include("./control/test_controlsystemgrowth.jl")
include("./control/test_controlsystemlinearized.jl")
include("./symbolic/test_automaton.jl")
include("./symbolic/test_symbolicmodel.jl")
include("./control/test_fromcontrolsystemgrowth.jl")
include("./control/test_fromcontrolsystemlinearized.jl")
include("./control/test_controller.jl")
include("./control/test_controllerreach.jl")
include("./control/test_controllersafe.jl")
