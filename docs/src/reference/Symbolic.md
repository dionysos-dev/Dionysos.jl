# Symbolic 

This folder contains the data structures needed to encode the different abstractions.

```@docs
Dionysos.Symbolic.AbstractSymbolicModel
Dionysos.Symbolic.SymbolicModel
Dionysos.Symbolic.GridBasedSymbolicModel
```

## Automaton
```@docs
Dionysos.Symbolic.compute_post!
```

```@docs
Dionysos.Symbolic.SymbolicModelList
Dionysos.Symbolic.LocalGridBasedSymbolicModel
```

```@docs
Dionysos.Symbolic.determinize_symbolic_model
```

```@docs
Dionysos.Symbolic.AbstractExecutionBackend
Dionysos.Symbolic.SequentialBackend
Dionysos.Symbolic.ThreadedBackend
Dionysos.Symbolic.JuliaDistributedBackend
Dionysos.Symbolic.SlurmArrayBackend
```