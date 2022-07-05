

# Manual of Dionysos
 This manual is still in the making.
## General Overview

Dionysos tackles, in a composite setting, optimal control/reachability/safety problems for complex dynamical systems, possibly in a black-box fashion.
For the sake of generality, the main components of Dionysos structure are highlighted in the following list, and detailed in specific subsections.


1. $\textit{Domain}$: Specifies the type of states/sets the systems will evolve in. Examples: Discrete sets, Hyperrectangles, Ellispioid 
2. $\textit{System}$: Defines system dynamics (on a domain type), growth bound, etc.
3. $\textit{Problem}$: Collection of System +  Objective function/reachability/target/LTL specification + Source/initial condition + Cost of the transiction + Additional Constraints.
4. $\textit{Solver}$: Defines  (abstraction-based) control design technique. Input: Problem. Output: Controller. 
5. $\textit{Controller}$: Control policy solving the problem.

##  Domain

## System


## Problem

## Solver

## Controller