# How to run the code 

## Example 1 : Default Controller 
In this example, you can find 2 codes : 

*  [`default_controller.jl`](1.%20Default%20Controller/default_controller.jl): Run the default ZMP based controller in order to save the reference joints trajectories in a [csv file](1.%20Default%20Controller/results/walkingPattern_ref.csv).
*  [`simulation_environment.jl`](1.%20Default%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a default ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

The defaut controller use the [param.jl file](/../deps/param.jl) in order to generate the trajectory. 

## Example 2 : Optimisation Process 
In this example, you can find 2 codes : 

*  [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl): Run an optimisation process by simulating all the candidates for a given discrtetised input-space. 
*  [`optimisation_post_process.jl`](2.%20Optimisation%20process/optimisation_post_process.jl): Manual post-process of the data collected using the [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl). 

## Example 3 : Optimised Controller 
In this example, you can find 2 codes : 

*  [`optimised_controller.jl`](3.%20Optimised%20Controller/optimised_controller.jl): Run the optimised ZMP based controller in order to save the reference joints trajectories in a [csv file](3.%20Optimised%20Controller/results/walkingPattern_ref.csv).
*  [`simulation_environment.jl`](3.%20Optimised%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a optimised ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

Once we found the optimised parameters using [Example 2 : Optimisation Process](2.%20Optimisation%20process/), the controller will not use the [param.jl file][../deps/param.jl] (the code will use them as initial values). 

## Example 4 : Simulation Environment 
In this example, you can find 1 code : 

*  [`simulation_environment.jl`](4.%20Simulation%20Environment/simulation_environment.jl): Open a csv file and fetch it to the simulation environment in order to simulate the trajectories given by the csv file. 

## Example 5 : Closed-loop Simulation 
Soon... 