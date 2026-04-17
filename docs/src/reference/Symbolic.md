# Symbolic 

This folder contains the data structures needed to encode the different abstractions.

```@docs
Dionysos.Symbolic.SymbolicModel
Dionysos.Symbolic.GridBasedSymbolicModel
```

```@docs
Dionysos.Symbolic.SymbolicModelList
Dionysos.Symbolic.LocalGridBasedSymbolicModel
```

```@docs
Dionysos.Symbolic.determinize_symbolic_model
```

```@docs
Dionysos.Symbolic.init_abstraction_workers!
Dionysos.Symbolic._install_distributed_abstraction_data!
Dionysos.Symbolic._warmup_distributed_abstraction_worker!
Dionysos.Symbolic.partition_source_state_ids
Dionysos.Symbolic._run_local_partition_ids
Dionysos.Symbolic.assign_states_to_workers
Dionysos.Symbolic._make_local_symmodel_from_ids
Dionysos.Symbolic._clear_distributed_abstraction_data!
Dionysos.Symbolic.clear_abstraction_workers!
```