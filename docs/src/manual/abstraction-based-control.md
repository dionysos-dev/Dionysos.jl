# Abstraction-based control

Given a mathematical description of the system dynamics and the specifications describing the desired closed-loop behavior of the system, 
**abstraction-based control** techniques involve synthesizing a **correct-by-construction** controller through a **systematic** three-step procedure. 
First, both the original system and the specifications are transposed into an abstract domain, resulting in an abstract system and corresponding abstract specifications. 
This step generally involves a complete discretization of the state and input spaces, typically with uniform hyperrectangles.
Next, an abstract controller is synthesized to solve this abstract control problem. Finally, the third step, referred to as the **refinement procedure**, involves deducing a controller for the original control problem from the abstract controller. The value of this approach lies in the substitution of the original system (often an infinite system) with a finite system, which enables it to leverage powerful control tools in the domain of symbolic control. This three steps procedure is illustrated on the following figure.

<figure>
    <center><img
    src="https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/assets/abstraction.png?raw=true"
    alt="Abstraction-based control."></center>
    <center><figcaption>Abstraction-based control</figcaption></center>
</figure>

Although this approach offers a safety-critical framework, it suffers from the curse of dimensionality due to the exponential growth of the number of states with respect to the dimension.
In order to render these techniques practical, it is necessary to construct **smart abstractions**, i.e., they differ from classical techniques in that the partitioning is designed smartly, using optimization-based design techniques, and computed iteratively, unlike the classical approach which uses an a priori defined approach, sub-optimal and subject to the curse of dimensionality.
To this end, we introduce solvers called **lazy solvers** (i.e. postponing heavier numerical operations) that co-design the abstraction and the controller to reduce the computed part of the abstraction.