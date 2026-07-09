# Dionysos
<picture>
  <source srcset="assets/images/logo-dark.png" media="(prefers-color-scheme: dark)">
  <img src="assets/images/logo.png"  height="240">
</picture>

| **Documentation & Paper** | **Build Status** | **Quality** |
|:-----------------:|:----------------:|:------------:|
| [![DOI][paper-img]][paper-url] [![][docs-latest-img]][docs-latest-url] [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] [![Codecov][codecov-img]][codecov-url] | [![Aqua QA][aqua-img]][aqua-url] |


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

[aqua-img]: https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

## Overview

This branch contains the implementation of the interactive natural-language interface presented in the paper **"End-to-End Abstraction-Based Control with LLM-Enhanced NL-to-LTL Translation"**. The interface translates natural-language requirements into Linear Temporal Logic (LTL), validates the generated specification with the user, synthesizes a controller, and visualizes the resulting robot motion.

## Installation

Install Julia and follow the package installation instructions described in the Julia documentation.

In addition, install **GLMakie** and **Genie**:

```julia
import Pkg
Pkg.add("GLMakie")
Pkg.add("Genie")
```

To use the LLM, create a `.env` file in the root of the repository and add your API key for the corresponding provider (e.g., Anthropic).

## Running the Interface

The main entry point is:

```text
InteractiveWithFeedBack.jl
```

Run this file from Julia. Once the server starts, open your web browser and navigate to:

```text
http://127.0.0.1:8000/
```

You can then enter natural-language requirements, validate the generated LTL specification, and synthesize and execute the corresponding controller through the web interface.
