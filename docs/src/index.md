```@raw html
<img src="assets/logo.png" height="240">
```
# Introduction

Welcome to the documentation for Dionysos!

## What is Dionysos ?

[Dionysos](https://github.com/dionysos-dev/Dionysos.jl) is the software of the ERC project [Learning to control](https://perso.uclouvain.be/raphael.jungers/content/erc-consolidator-grant) (L2C) embedded in [Julia](https://julialang.org/). In view of the Cyber-Physical Revolution, the only sensible way of controlling these complex systems is often by discretizing the different variables, thus transforming the model into a simple combinatorial problem on a finite-state automaton, called an abstraction of this system. The goal of L2C is to transform this approach into an effective, scalable, cutting-edge technology that will address the CPS challenges and unlock their potential. This ambitious goal will be achieved by leveraging powerful tools from Mathematical Engineering.

## Current version

The current version is still in the making, and allows to solve problems such as reachability problems for hybrid systems.

## Longterm objectives
Rather than relying on closed-form analysis of a model of the dynamical system, Dionysos will learn the optimal control from data, whether harvested from the physical system or generated synthetically. It will rely on a novel methodology, combining the efficiency of several modern optimization/control-theoretic/machine-learning techniques with the theoretical power of the Abstraction approach. All the pieces of the architecture are chosen to foster black-box and data-driven analysis, thereby matching rising and unresolved challenges. Summarizing, the objectives are
* To develop a mathematical and algorithmic framework for efficient Abstraction of Cyber-Physical Systems thriving on recent technologies in Optimization and Control;
* To leverage this framework in situations where the system is described by data, rather than a classical model.

## Structure of the documentation 

The documentation is organised as follows.

* The **Examples** section contains a few examples of solving problems with Dionysos. Start with [Getting started](https://dionysos-dev.github.io/Dionysos.jl/dev/generated/Getting%20Started/) if you want to get familiar with Dionysos.
* The **Manual** section contains all the useful information to use Dionysos as a user.
* The **API Reference** sections contains all the functions that we currently use in Dionysos. 
* The **Developer Docs** section is dedicated to the contributors to Dionysos developement. 

 
## Need help?

If you need help, open an [issue](https://github.com/dionysos-dev/Dionysos.jl/issues) on Github.

## ERC sponsor 

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme under grant agreement No 864017 - L2C.

```@raw html
<img class="display-light-only" src="assets/logo_erc_white.jpg" alt="ERC logo"/>
<img class="display-dark-only" src="assets/logo_erc_black.jpg" alt="ERC logo"/>
```