# How to install and run the 7380xing/Dionysos.jl project 

* Install Julia and Git Bash
* Open Git Bash and go to your desired directory `~/`
* Use the following command `git clone https://github.com/7380Xing/Dionysos.jl.git` to clone the repository 
```
$ git clone https://github.com/7380Xing/Dionysos.jl.git
Cloning into 'Dionysos.jl'...
remote: Enumerating objects: 3664, done.
remote: Counting objects: 100% (1662/1662), done.
remote: Compressing objects: 100% (870/870), done.
remote: Total 3664 (delta 981), reused 1317 (delta 777), pack-reused 2002
Receiving objects: 100% (3664/3664), 37.48 MiB | 426.00 KiB/s, done.
Resolving deltas: 100% (2219/2219), done.
Updating files: 100% (212/212), done.
```
* Then, go to the clonned directory 
```
$ cd Dionysos.jl/
```
## Run on Visual Studio Code 
* Run the code `BipedRobot//ZMPBipedRobot//examples//1. Default Controller//default_controller.jl` in order to launch Julia REPL
* This will lead to an error but we need to activate the `ZMPBipedRobot` environment in this new Julia REPL terminal 
```
julia> ]
(@v1.9) pkg> activate BipedRobot\\ZMPBipedRobot\\
  Activating project at `E:\Design and control of a biepiedal walking robot\temp\7380xing\Dionysos.jl\BipedRobot\ZMPBipedRobot`
(ZMPBipedRobot) pkg> 
```
* Then add all necessary pacakges in the new environment using `instantiate` command
```
(ZMPBipedRobot) pkg> instantiate
```
* Wait until all packages is installed
* Run the code `BipedRobot//ZMPBipedRobot//test//runtest.jl` as a test 

* Run the code in the following folder for more examples

## Example 1 : Default Controller 
In this example, you can find 2 codes : 

*  [`default_controller.jl`](1.%20Default%20Controller/default_controller.jl): Run the default ZMP based controller in order to save the reference joints trajectories in a [csv file](../docs/1.%20Default%20Controller/walkingPattern_ref.csv).
*  [`simulation_environment.jl`](1.%20Default%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a default ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

The defaut controller use the [param.jl file](/../deps/param.jl) in order to generate the trajectory. 

## Example 2 : Optimisation Process 
In this example, you can find 2 codes : 

*  [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl): Run an optimisation process by simulating all the candidates for a given discrtetised input-space. 
*  [`optimisation_post_process.jl`](2.%20Optimisation%20process/optimisation_post_process.jl): Manual post-process of the data collected using the [`optimisation_process.jl`](2.%20Optimisation%20process/optimisation_process.jl).

Note : The code will take a lot of time (> 1 hour) because the "optimisation process" just compute for every candidate their objective values. 

## Example 3 : Optimised Controller 
In this example, you can find 2 codes : 

*  [`optimised_controller.jl`](3.%20Optimised%20Controller/optimised_controller.jl): Run the optimised ZMP based controller in order to save the reference joints trajectories in a [csv file](/../docs/3.%20Optimised%20Controller/walkingPattern_ref_slow.csv). 
*  [`simulation_environment.jl`](3.%20Optimised%20Controller/simulation_controller.jl): Generate the reference joints trajectories from a optimised ZMP based controller, in order to fetch it to the simulation environment and run the simulation. 

Once we found the optimised parameters using [Example 2 : Optimisation Process](2.%20Optimisation%20process/) (by using the variable `best_key`), the controller will not use the [param.jl file](../deps/param.jl) (the code will use them as initial values). 

## Example 4 : Simulation Environment 
In this example, you can find 1 code : 

*  [`simulation_environment.jl`](4.%20Simulation%20Environment/simulation_environment.jl): Open a csv file and fetch it to the simulation environment in order to simulate the trajectories given by the csv file. 

To run this code, you need first generate the CSV file from  [`default_controller.jl`](1.%20Default%20Controller/default_controller.jl) or  [`optimised_controller.jl`](3.%20Optimised%20Controller/optimised_controller.jl) and change the folder path accordingly. 

Note : This code is only used for a specific type of robot which require flat foot. Other type of URDF file which doesn't follow the same structure as the example may not work on this code. 
