# Symbolic 

This folder contains the data structures needed to encode the different abstractions:
the finite automaton (transition graph), the symbolic models built on top of a mapping,
the parallel execution backends that populate them, and the timed hybrid symbolic model.

## Symbolic models

```@docs
Dionysos.Symbolic.AbstractSymbolicModel
Dionysos.Symbolic.SymbolicModel
Dionysos.Symbolic.GridBasedSymbolicModel
Dionysos.Symbolic.SymbolicModelList
Dionysos.Symbolic.LocalGridBasedSymbolicModel
Dionysos.Symbolic.determinize_symbolic_model
```

## Automaton

The transition graph. Concrete implementations trade memory for `pre`/`post` speed.

```@docs
Dionysos.Symbolic.AbstractAutomatonList
Dionysos.Symbolic.SortedAutomatonList
Dionysos.Symbolic.IndexedAutomatonList
Dionysos.Symbolic.FastIndexedAutomatonList
Dionysos.Symbolic.compute_post!
Dionysos.Symbolic.finalize!
```

## Transition metadata

```@docs
Dionysos.Symbolic.TransitionKey
```

## Execution backends

How the transition relation is computed from a concrete system.

```@docs
Dionysos.Symbolic.AbstractExecutionBackend
Dionysos.Symbolic.SequentialBackend
Dionysos.Symbolic.ThreadedBackend
Dionysos.Symbolic.JuliaDistributedBackend
Dionysos.Symbolic.SlurmArrayBackend
```

## Timed hybrid symbolic model

Abstraction of a timed hybrid system: per-mode spatial + time abstractions flattened into
one automaton, with inputs unified through a global input map. The optimizer-driven builders
live in `Optim/hybrid_systems/HybridSystemAbstraction`; the data structures live here.

```@docs
Dionysos.Symbolic.TimedHybridSymbolicModel
Dionysos.Symbolic.TimeSymbolicModel
Dionysos.Symbolic.GlobalInputMap
```
