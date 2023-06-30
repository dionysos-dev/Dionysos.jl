# Two following tests temporarily removed, see https://github.com/dionysos-dev/Dionysos.jl/issues/209
# include("./utils/data_structures/BDD/test_BDD.jl")
# include("./utils/data_structures/BDD/test_inttupleset.jl")

include("./utils/data_structures/test_queue.jl")
include("./utils/data_structures/test_digraph.jl")
include("./utils/data_structures/tree.jl")
include("./utils/search/test_search.jl")
include("./utils/rectangle.jl")
include("./utils/ellipsoid.jl")
include("./utils/scalar_functions.jl")

include("./domain/test_griddomain.jl")
include("./domain/test_general_domain.jl")
include("./domain/test_nested_domain.jl")

include("./symbolic/test_automaton.jl")
include("./symbolic/test_proba_automaton.jl")
include("./symbolic/test_symbolicmodel.jl")
include("./symbolic/test_lazy_symbolic.jl")
include("./symbolic/test_ellipsoidal_transitions.jl")

include("./problem/test_problems.jl")

include("./control/test_controller.jl")
include("./control/test_controllerreach.jl")
include("./control/test_controllersafe.jl")
include("./control/test_fromcontrolsystemgrowth.jl")
include("./control/test_fromcontrolsystemlinearized.jl")

include("./optim/test_SCOTS_safety.jl")
include("./optim/test_SCOTS_reachability.jl")
include("./optim/test_lazy_abstraction.jl")
include("./optim/test_ellipsoids_abstraction.jl")
include("./optim/test_lazy_ellipsoids_abstraction.jl")
include("./optim/test_hierarchical_abstraction.jl")

include("./mapping/test_mapping_continuous.jl")

include("./system/test_controlsystem.jl")
include("./system/test_controller.jl")

include("./examples/test_gol_lazar_belta.jl")
