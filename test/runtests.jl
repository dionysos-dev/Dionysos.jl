include("./utils/data_structures/BDD/test_BDD.jl")
include("./utils/data_structures/BDD/test_inttupleset.jl")
include("./utils/data_structures/test_queue.jl")
include("./utils/search/test_search.jl")
include("./utils/rectangle.jl")

include("./domain/test_griddomain.jl")
include("./domain/test_general_domain.jl")
include("./domain/test_nested_domain.jl")
include("./domain/test_general_domain.jl")


include("./symbolic/test_automaton.jl")
include("./symbolic/test_symbolicmodel.jl")
include("./symbolic/test_lazy_symbolic.jl")

include("./control/test_fromcontrolsystemgrowth.jl")
include("./control/test_fromcontrolsystemlinearized.jl")
include("./control/test_controller.jl")
include("./control/test_controllerreach.jl")
include("./control/test_controllersafe.jl")
include("./control/test_lazy_periodicity.jl")


include("./examples/test_state_feedback_trans.jl")
include("./examples/test_gol_lazar_belta.jl")
