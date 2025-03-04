# Dionysos
<picture>
  <source srcset="assets/images/logo-dark.png" media="(prefers-color-scheme: dark)">
  <img src="assets/images/logo.png"  height="240">
</picture>

| **Documentation and paper** | **Build Status** |
|:-----------------:|:----------------:|
| [![DOI][paper-img]][paper-url] [![][docs-latest-img]][docs-latest-url]  [![][docs-stable-img]][docs-stable-url]        | [![Build Status][build-img]][build-url]  [![Codecov branch][codecov-img]][codecov-url]      |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://dionysos-dev.github.io/Dionysos.jl/stable
[docs-latest-url]: https://dionysos-dev.github.io/Dionysos.jl/dev
[paper-img]: https://proceedings.juliacon.org/papers/10.21105/jcon.00160/status.svg
[paper-url]: https://doi.org/10.21105/jcon.00160

[build-img]: https://github.com/dionysos-dev/Dionysos.jl/actions/workflows/ci.yml/badge.svg?branch=master
[build-url]: https://github.com/dionysos-dev/Dionysos.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/github/dionysos-dev/Dionysos.jl/coverage.svg
[codecov-url]: https://app.codecov.io/github/dionysos-dev/Dionysos.jl

## Overview
Dionysos is the software of the ERC project Learning to control (L2C). In view of the Cyber-Physical Revolution, the only sensible way of controlling these complex systems is often by discretizing the different variables, thus transforming the model into a simple combinatorial problem on a finite-state automaton, called an abstraction of this system. The goal of L2C is to transform this approach into an effective, scalable, cutting-edge technology that will address the CPS challenges and unlock their potential. This ambitious goal will be achieved by leveraging powerful tools from Mathematical Engineering.

## Current version
The current version is still in the making, and allows to solve problems such as reachability problems for hybrid systems. See the [Docs](https://dionysos-dev.github.io/Dionysos.jl/dev/) for further information.

## Longterm objectives
Rather than relying on closed-form analysis of a model of the dynamical system, Dionysos will learn the optimal control from data, whether harvested from the physical system or generated synthetically. It will rely on a novel methodology, combining the efficiency of several modern optimization/control-theoretic/machine-learning techniques with the theoretical power of the Abstraction approach. All the pieces of the architecture are chosen to foster black-box and data-driven analysis, thereby matching rising and unresolved challenges. Summarizing, the objectives are
* To develop a mathematical and algorithmic framework for efficient Abstraction of Cyber-Physical Systems thriving on recent technologies in Optimization and Control;
* To leverage this framework in situations where the system is described by data, rather than a classical model;

## Installation

Download Julia, and follow the instructions described [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-packages).
