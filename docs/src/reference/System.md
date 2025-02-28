# System

This folder contains different ways to define systems, for instance to encode a controller.

## Control system

Each control system should be implemented as a ControlSystem.

```@docs
Dionysos.System.ControlSystem
```

So far, we have implemented a few examples of control systems: 

```@docs
Dionysos.System.SimpleSystem
Dionysos.System.ControlSystemGrowth
Dionysos.System.ControlSystemLinearized
Dionysos.System.EllipsoidalAffineApproximatedSystem
```

## Controller 
So far, the abstraction-based methods that we use define either piecewise-constant or piecewise-affine controllers.

```@docs
Dionysos.System.ConstantController
Dionysos.System.AffineController
```

## Trajectories 
```@docs
Dionysos.System.DiscreteTrajectory
Dionysos.System.ContinuousTrajectory
Dionysos.System.HybridTrajectory
Dionysos.System.Trajectory
Dionysos.System.Control_trajectory
Dionysos.System.Cost_control_trajectory
```

## Approximation 
```@docs
Dionysos.System.get_under_approximation_map
Dionysos.System.get_over_approximation_map
```