# Dionysos

<img src="logo.png" height="240">

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://dionysos-dev.github.io/Dionysos.jl/stable
[docs-latest-url]: https://dionysos-dev.github.io/Dionysos.jl/dev

[build-img]: https://github.com/dionysos-dev/Dionysos.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/dionysos-dev/Dionysos.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/dionysos-dev/Dionysos.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/dionysos-dev/Dionysos.jl?branch=master

## Overview
Dionysos is the software of the ERC project Learning to control (L2C). In view of the Cyber-Physical Revolution, the only sensible way of controlling these complex systems is often by discretizing the different variables, thus transforming the model into a simple combinatorial problem on a finite-state automaton, called an abstraction of this system. The goal of L2C is to transform this approach into an effective, scalable, cutting-edge technology that will address the CPS challenges and unlock their potential. This ambitious goal will be achieved by leveraging powerful tools from Mathematical Engineering.

## Current version

The current version is still in the making, and allows to solve problems such as reachability problems for hybrid systems. See the [Examples](https://github.com/dionysos-dev/Dionysos.jl/tree/master/examples) for further information.

## Longterm objectives
Rather than relying on closed-form analysis of a model of the dynamical system, Dionysos will learn the optimal control from data, whether harvested from the physical system or generated synthetically. It will rely on a novel methodology, combining the efficiency of several modern optimization/control-theoretic/machine-learning techniques with the theoretical power of the Abstraction approach. All the pieces of the architecture are chosen to foster black-box and data-driven analysis, thereby matching rising and unresolved challenges. Summarizing, the objectives are
* To develop a mathematical and algorithmic framework for efficient Abstraction of Cyber-Physical Systems thriving on recent technologies in Optimization and Control;
* To leverage this framework in situations where the system is described by data, rather than a classical model;

## Installation

Download Julia, and follow the instructions described [here](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-unregistered-packages).
