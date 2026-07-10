# Utils 

This folder contains all the auxiliary functions needed.

## Functions 

```@docs
Dionysos.Utils.ScalarFunction
Dionysos.Utils.ScalarControlFunction
Dionysos.Utils.QuadraticFunction
Dionysos.Utils.QuadraticStateControlFunction 
Dionysos.Utils.PolyhedralFunction
Dionysos.Utils.BlackBoxFunction
Dionysos.Utils.BlackBoxControlFunction
```

## Search

```@docs
Dionysos.Utils.RRT
Dionysos.Utils.Tree
Dionysos.Utils.NodeT
Dionysos.Utils.add_node!
```

## Geometric shapes

```@docs
Dionysos.Utils.Box
Dionysos.Utils.box
Dionysos.Utils.get_quadratic_form
Dionysos.Utils.get_sublevel_set
Dionysos.Utils.get_length_semiaxis
Dionysos.Utils.set_in_period
Dionysos.Utils.get_min_bounding_box
Dionysos.Utils.project_set
Dionysos.Utils.sample
Dionysos.Utils.samples
```

## Set algebra

```@docs
Dionysos.Utils.set_union
Dionysos.Utils.set_minus
Dionysos.Utils.minus_included
Dionysos.Utils.minus_hole
Dionysos.Utils.empty_region
Dionysos.Utils.is_included
Dionysos.Utils.is_disjoint
```

## Discretization helpers

```@docs
Dionysos.Utils.invert_incl_mode
Dionysos.Utils.wrap_value
```

## Scalar optimization

```@docs
Dionysos.Utils.golden_section_search
Dionysos.Utils.newton_method
Dionysos.Utils.derivative_bisection
```


## Path-Complete Framework

```@docs
Dionysos.Utils.PathCompleteFramework.LabDigraph
Dionysos.Utils.PathCompleteFramework.PCLF
Dionysos.Utils.PathCompleteFramework.compute_quadratic_pieces_pclf
Dionysos.Utils.PathCompleteFramework.compute_symmetric_2n_faces_polyhedral_pieces_pclf
Dionysos.Utils.PathCompleteFramework.compute_polyhedral_pieces_pclf
```
