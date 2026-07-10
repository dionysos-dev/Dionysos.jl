# Utils 

This folder contains all the auxiliary functions needed.

## Functions 

```@docs
Dionysos.Utils.QuadraticStateControlFunction 
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
Dionysos.Utils.HyperRectangle
Dionysos.Utils.DeformedRectangle
Dionysos.Utils.set_in_period
Dionysos.Utils.get_min_bounding_box
```

## Set algebra

```@docs
Dionysos.Utils.set_union
Dionysos.Utils.set_minus
Dionysos.Utils.minus_included
Dionysos.Utils.minus_hole
Dionysos.Utils.empty_region
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
