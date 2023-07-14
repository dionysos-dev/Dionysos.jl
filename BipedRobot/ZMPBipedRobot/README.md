# Design and Control of a Biped Robot 

## Implementation of an ZMP based controller for a simple 8 D.o.F robot 

This project was part of another project called [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) in order to benchmark the biped robot. This project is associated to the following  [thesis](https://eur03.safelinks.protection.outlook.com/?url=https%3A%2F%2Fdial.uclouvain.be%2Fdownloader%2Fdownloader_thesis.php%3Fpid%3Dthesis%3A40693%26datastream%3DPDF_01%26key%3D8b2cc138cd5db26d48602e804a9a548a&data=05%7C01%7Cxing.chen%40student.uclouvain.be%7C6c13df580a2143b3df6b08db6756670c%7C7ab090d4fa2e4ecfbc7c4127b4d582ec%7C0%7C0%7C638217393017003936%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=Z%2B0DtvalPfAVmtWYsjAwp0SMeM%2B058%2BAoZ2wLHBK7h0%3D&reserved=0
).
Every part of this project is sample code which shows how to do the following : 

* Generate the joint trajectories for a biped robot based on the ZMP stability criteria on a csv file 
* Read an URDF file 
* Create a robot visualiser 
* Simulate the biped robot on a well-constrolled environment 

## Controlle r's Structure 
![Structure of the ZMP controller](ZMP_based_controller.png)

The figure shows the intended structure of the implemented controller in open loop. The controller separated into 2 mains blocks : 
* ZMP based controller : 
    * Foot Planner, evaluates the landing position of the right and left foot.
    * Swing Foot Trajectory, determines the 3D position of the swing foot.
    * ZMP Trajectory Generator, defines the reference ZMP trajectory to remain within the support polygon.
    * CoM Trajectory Generator, computes the CoM trajectory and the hip body position based on the reference ZMP.
    * Inverse Kinematics, converts the local foot position in workspace coordinates intojoint space coordinates. 
* Simulation Environment : 
    * Position Control, a classical PID controller with dynamic compensation.
    * Robot, a virtual robot in a virtual environment 
    
![Result of the ZMP controller](simres_opti_robot.gif)

## How to run this project 

This project has many examples, see [Examples](examples/) for further information.

## Actual Version 

This version of the project does not support a closed-loop system. In a short term, a closed-loop form will be develloped in a disturbed environment. 