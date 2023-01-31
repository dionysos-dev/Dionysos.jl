

# **Manual of Dionysos**
 This manual provides a brief overview of the framework, problems and solving strategies tackled by Dionysos. 



 # High-Level Description
Dionysos tackles, in a composite setting, optimal control/reachability/safety problems for complex dynamical systems, possibly in a black-box and data-driven fashion.
In what follows, we provide a description of the main steps Dionyosos follows in order to encode and possibly solve complex control problems on hybrid systems. 
## Initialization

This initial step concerns all the aspects related to the initialization of the main algorithm.

 It should be fed with 
 1.   A proper representation/definition of the domain state and relative constraints,
 2.   Structure and system dynamics (in closed form or data-driven fashion), 
 3.   Control/optimization/reachability/safety tasks or objectives, from now on called ``the problem''. 

In particular, in this initialization step of the algorithm, a suitable system representation should be chosen, either in terms of data (in the case of data-driven design) or differential/difference equations or hybrid systems, with the possibility of having a combination of both paradigms.

The problem formulation formalizes the tasks we want to solve, that is, given a suitable representation of the system dynamics (possible with limited information), the problem clearly states the collection of specifications that should be met by the trajectories of the system obtained applying the control policy provided by the algorithm. These specifications may include aspects related to closed-loop system performance (e.g. safety/reachability/stability ), cost of implementation (via objective functions), computational complexity, etc.

## Main algorithm
The main algorithm takes the system representation and the problem formulation as inputs and returns a suitable solution, in terms of either a certificate of feasibility (and the associated control policy) or a KO/unfeasability message. 

It has been designed in a modular flavor, in order to ensure that even important changes to one module would not harm the overall structure. The major advantage of this approach lies in the possibility to integrate new or pre-existing solvers, to form a global solver capable of tackling complex control problems.

In order to achieve this modularity, we have to identify generic modules dedicated to solve some particular sub-tasks. Then, we need an automatic orchestrator, that, depending on the particular system and problem, under consideration, chooses which modules activate and in which order.

### **Orchestrator**
Depending on the system specifications, a suitable combination of different strategies may be chosen to reach a solution. The orchestrator is a meta-solver whose automatically (but, possibly, under a partial guide of the user) chooses the order of activation of different solvers. 
The user can thus choose to specify a pre-decided control strategy (i.e. a specific sub-module) or to let the orchestrator decide and test different strategy, or a combination of the two.


The orchestrator supervises the sub-tasks being run, both in a global and local level, and taking decisions on which next steps should be taken depending on the current state of the program in respect to the specifications. 
 
 ### **Global solvers: Abstractions**
The main global technique employed by Dionysos is an *Abstraction based approach*: The input system *S* and the corresponding control problem *P* are ``translated'' into corresponding finite representations *S'* and *P'*, together with a map *F* (a simulation relation), formally relating (solutions of) the 2 systems and corresponding control tasks/cost. 

Possible kind of abstractions are described in what follows:
- Hyper-graphs based abstraction: the domain is split into a finite number of sub-cells, and the corresponding finite-state system provides a non-deterministic abstraction, formalizing a worst-case analysis approach;
- Markov Decision Process (MDP): each abstract transition inherits a probability, possibly arising from the data/statistical-based definition of the system;
-  Interval Markov Decision process: taking .

This finite global representation is then used to infer global information on the original control problem. Global techniques on the abstracted systems are, for examples:

- Shortest path techniques;
-   Value/Lyapunov function iterations;
-  Steady-state vector, entropy of MDP,
-  Model predictive control;
- Path-complete techniques;

These global techniques are considered as ''modules'' by the orchestrator and can be used sequentially or even in parallel.

 ### **Local solvers/iterations**

Based on global information obtained from intermediate steps/abstractions, the orchestrator/user can decide to take local actions to improve the quality of the solution of the original control problem.

For example, possible local approaches could impose:
-  To refine the abstraction by cutting cells, merging cells, possibly only locally;
-  To compute the abstraction in parts of the state space which were not explored in the previous global step;
-  To compute/improve local controllers in a subset of the (abstract) state space to reduce the cost or the non-determinism level of the current abstraction;
-  To collect additional data sampling the dynamics in some parts of the state space (in a data-driven approach).






# General Overview of the Structures
<span style="color:blue">This section will possibly change, with new possible objects/strucutres text</span>

The main components of Dionysos structure are highlighted in the following list, and detailed in specific subsections.



1. $\textit{Domain}$: Specifies the type of states/sets the systems will evolve in. Examples: Discrete sets, Hyperrectangles, Ellispioids 
2. $\textit{System}$: Defines system dynamics (on a domain type), growth bound, etc. General enough to tackle hybrid systems (joint continuous and discrete time behaviour)
3. $\textit{Problem}$: Collection of System +  Objective function/reachability/target/LTL specification + Source/initial condition + Cost of the transiction + Additional Constraints.
4. $\textit{Solver}$: Defines  (abstraction-based) control design technique. Input: Problem. Output: Controller. 
5. $\textit{Controller}$: Control policy solving the problem.

##  Domain

## System


## Problem

## Solver

## Controller