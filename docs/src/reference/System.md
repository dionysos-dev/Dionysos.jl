# System

This folder contains different ways to define systems, for instance to encode a controller.

## Control system

Each control system should be implemented as a ControlSystem .

```@docs
Dionysos.System.ControlSystem
```

So far, we have implemented a few examples of control systems : 

```@docs
Dionysos.System.SimpleSystem
Dionysos.System.ControlSystemGrowth
Dionysos.System.ControlSystemLinearized
Dionysos.System.EllipsoidalAffineApproximatedSystem
```

## Controller 
So far, the abstraction-based methods that we use define either piece-wise constant or piecewise-affine controllers.

```@docs
Dionysos.System.ConstantController
Dionysos.System.AffineController
```

## Trajectories 
```@docs
Dionysos.Control.DiscreteTrajectory
Dionysos.Control.ContinuousTrajectory
```