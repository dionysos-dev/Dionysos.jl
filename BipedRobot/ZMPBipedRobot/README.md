# Design and Control of a Biped Robot 

## Implementation of an ZMP based controller for a simple 8 D.o.F robot 

This project is part of another project called [Dionysos](https://github.com/dionysos-dev/Dionysos.jl) in order to benchmark the biped robot. This project is associated to the following  [thesis](https://dial.uclouvain.be/downloader/downloader_thesis.php?pid=thesis:40693&datastream=PDF_01&key=8b2cc138cd5db26d48602e804a9a548a
).
Every part of this project is sample code which shows how to do the following : 

* Generate the joint trajectories for a biped robot based on the ZMP stability criteria on a csv file 
* Read an URDF file 
* Create a robot visualiser 
* Simulate the biped robot on a well-controlled environment 

![Robot 2D model](https://github.com/7380Xing/Dionysos.jl/assets/99494151/46b26ba4-53af-4dd0-936a-76c8f2c6e123)


## Controller Structure 
![Structure of the ZMP controller](https://github.com/7380Xing/Dionysos.jl/assets/99494151/4ddd3f21-071b-4264-a485-be3cac7fc3c5)
ark.png)

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
    
![Result of the ZMP controller](https://github.com/7380Xing/Dionysos.jl/assets/99494151/1112c75a-d8aa-47c2-9f44-c9a1254466fb)

## Main References 
| Block | Reference(s) |
|-------|--------------|
| Foot Planner | R. Khusainov, A. Sagitov, A. Klimchik, and E. Magid. “Arbitrary Trajectory Foot Planner for Bipedal Walking:” in: Proceedings of the 14th International Conference on Informatics in Control, Automation and Robotics. 14th International Conference on Informatics in Control, Automation and Robotics. Madrid, Spain: SCITEPRESS - Science and Technology Publications, 2017, pp. 417–424. isbn: 978-989-758-263-9 978-989-758-264-6. doi: 10.5220/ 006442504170424.|
|Swing Foot | R. Khusainov, A. Sagitov, A. Klimchik, and E. Magid. “Arbitrary Trajectory Foot Planner for Bipedal Walking:” in: Proceedings of the 14th International Conference on Informatics in Control, Automation and Robotics. 14th International Conference on Informatics in Control, Automation and Robotics. Madrid, Spain: SCITEPRESS - Science and Technology Publications, 2017, pp. 417–424. isbn: 978-989-758-263-9 978-989-758-264-6. doi: 10.5220/ 006442504170424.|
| CoM Trajectory Generator | S. Kajita, F. Kanehiro, K. Kaneko, K. Fujiwara, K. Harada, K. Yokoi, and H. Hirukawa. “Biped walking pattern generation by using preview control of zero-moment point”. In: 2003 IEEE International Conference on Robotics and Automation (Cat. No.03CH37422). IEEE International Conference on Robotics and Automation. IEEE ICRA 2003. Taipei, Taiwan: IEEE, 2003, pp. 1620–1626. isbn: 978-0-7803-7736-3. doi: 10.1109/ROBOT.2003.1241826.|
| Preview Control  | T. Katayama, T. Ohki, T. Inoue, and T. Kato. “Design of an optimal controller for a discrete-time system subject to previewable demand”. In: International Journal of Control 41.3 (Mar. 1985), pp. 677–699. issn: 0020-7179, 1366-5820. doi: 10 . 1080 / 0020718508961156.|
## How to run this project 

This project has many examples, see [Examples](examples/) for further information.

## Actual Version 

This version of the project does not support a closed-loop system. In a short term, a closed-loop form will be developed to handle disturbed environment. 
