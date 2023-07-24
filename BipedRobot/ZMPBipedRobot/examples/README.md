# How to run the code 

* Install Julia on your computer
* Activate the ZMPBipedRobot project by using the following command in the Julia REPL 

```
julia> ] 
(@v1.9) pkg> activate BipedRobot\\ZMPBipedRobot\\
    Activating project at `YourPath\Dionysos.jl\BipedRobot\ZMPBipedRobot`
(ZMPBipedRobot) pkg> 
```
* Run the code in the following folder 

## Example 1 : Default Controller 
In this example, you can find 2 codes : 

*  [`default_controller.jl`](1.%20Default%20Controller/default_controller.jl): Run the default ZMP based controller in order to save the reference joints trajectories in a [csv file](../docs/1.%20Default%20Controller/walkingPattern_ref.csv).
*  [`simulation_environment.jl`](1.%20Default%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a default ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

The defaut controller use the [param.jl file](/../deps/param.jl) in order to generate the trajectory. 

## Example 2 : Optimisation Process 
In this example, you can find 2 codes : 

*  [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl): Run an optimisation process by simulating all the candidates for a given discrtetised input-space. 
*  [`optimisation_post_process.jl`](2.%20Optimisation%20process/optimisation_post_process.jl): Manual post-process of the data collected using the [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl). 

## Example 3 : Optimised Controller 
In this example, you can find 2 codes : 

*  [`optimised_controller.jl`](3.%20Optimised%20Controller/optimised_controller.jl): Run the optimised ZMP based controller in order to save the reference joints trajectories in a [csv file](/../docs/3.%20Optimised%20Controller/walkingPattern_ref_slow.csv).
*  [`simulation_environment.jl`](3.%20Optimised%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a optimised ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

Once we found the optimised parameters using [Example 2 : Optimisation Process](2.%20Optimisation%20process/), the controller will not use the [param.jl file](../deps/param.jl) (the code will use them as initial values). 

## Example 4 : Simulation Environment 
In this example, you can find 1 code : 

*  [`simulation_environment.jl`](4.%20Simulation%20Environment/simulation_environment.jl): Open a csv file and fetch it to the simulation environment in order to simulate the trajectories given by the csv file. 

Note : This code is only used for a specific type of robot which require flat foot. Other type of URDF file which doesn't follow the same structure as the example may not work on this code. 

## Example 5 : Closed-loop Simulation 
Soon... 

# Example 6 : Generic Simulation Environment 

In this example, you can find 1 code : 

*  [`simulation_environment.jl`](6.%20Generic%20Simulation%20Environment/simulation_environment.jl): Open a csv file and an URDF file and fetch it to the simulation environment in order to simulate the trajectories given by the csv file. 

This code also need to define a controller which is done by defining the function `controller!(τ, t, state)` where  `τ` is the previous output torque at time `k - 1`, `k` is the actual sample, `t` is the actual time at `k` and `state` is the current state of the `mechanism` defined by `RigidBodyDynamics.jl`. This controller will return a torque which is directly ussed in the `RigidBodyDynamics.simulate` function.  
 