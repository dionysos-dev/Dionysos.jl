# Generic Simulation Environment 

In this example, you can find 1 code : 

*  [`simulation_environment.jl`](6.%20Generic%20Simulation%20Environment/simulation_environment.jl): Open a csv file and an URDF file and fetch it to the simulation environment in order to simulate the trajectories given by the csv file. 

This code also need to define a controller which is done by defining the function `controller!(τ, t, state)` where  `τ` is the previous output torque at time `k - 1`, `k` is the actual sample, `t` is the actual time at `k` and `state` is the current state of the `mechanism` defined by `RigidBodyDynamics.jl`. This controller will return a torque which is directly ussed in the `RigidBodyDynamics.simulate` function.  

The contact model needs to be difined as well. Here, we are using Hunt-Crossley model for the normal direction and Coulomb model for tangent direction.  

Other input are the location of the spherical contact points w.r.t the ankle frame. 

As an example, an basic 1-joint "bang-bang" controller is implemented which will oscillate between 2 values. The expected result is that the torque should change like a triangle wave form. 