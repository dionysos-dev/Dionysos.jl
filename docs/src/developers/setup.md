# Set up

This guide shows you what to do when you start developing for Dionysos.

## Installations

Start by installing Julia, VSCode and the Julia extension of VSCode as detailed [here](https://code.visualstudio.com/docs/languages/julia#_getting-started).
Now, install Git by following https://git-scm.com/book/en/v2/Getting-Started-Installing-Git.

## Launching the prompts

For every step, we show both how to do it from Visual Studio Code (**VSCode**) or from the **Julia REPL** or **Git bash**.

### Start the Julia REPL

To start a Julia REPL, type `Shift+Ctrl+P` and then `Julia: Start REPL`.
You should have installed Julia and the Julia VSCode extension as detailed in [Installations](@ref) for this to work.
You should see a prompt `julia>` appearing.
We always show the prompt you should see for every command as well as the output, don't copy-paste the prompt nor the output.

### Start Julia Pkg prompt

First [Start the Julia REPL](@ref).
Then, by pressing the `]` character you will see a `(@v1.8) pkg>` prompt appearing (if you are using the global environment). By pressing the backspace you will get back to the `julia>` prompt.
We always show the prompt you should see for every command as well as the output, don't copy-paste the prompt nor the output.

### Start Git bash

To start Git bash, click on the top menu of VSCode on "View" then "Terminal". On the bottom right, click on the down arrow at the right of the "+" and then on "Git bash" on the dropdown menu that appears.
You should have installed Julia and the Julia VSCode extension as detailed in [Installations](@ref) for this to work.
You should see a prompt `$` appearing.
We always show the prompt you should see for every command as well as the output, don't copy-paste the prompt nor the output.

## Cloning Dionysos

The purpose of this is to clone Dionysos at the location `~/.julia/dev/Dionysos` where `~` is your home folder.

### VSCode

Switch to Source Control by pressing `Ctrl+Shift+G` then on the three horizontal dots on the top right of the left pane then "clone"
then write `https://github.com/dionysos-dev/Dionysos.jl.git` and then select the folder `.julia/dev`.
Then rename the created folder `~/.julia/dev/Dionysos` into `~/.julia/dev/Dionysos` using your file manager.

### Julia REPL

See [Start the Julia REPL](@ref).

```julia
julia> using Pkg; Pkg.develop(url="https://github.com/dionysos-dev/Dionysos.jl.git")
```

### Julia Pkg prompt

See [Start Julia Pkg prompt](@ref).

```julia
(@v1.8) pkg> dev https://github.com/dionysos-dev/Dionysos.jl.git
```

## Install the Revise.jl and JuliaFormatter.jl packages

[Revise.jl](https://github.com/timholy/Revise.jl) reduces the need to restart your Julia REPL when you make changes in the source code. [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) allows to format your code following the rules stated in `.JuliaFormatter.toml`.

We install these two packages in the global environment so that it is available from all the environments.

### Julia REPL

See [Start the Julia REPL](@ref).

```julia
julia> using Pkg; Pkg.add("Revise"); Pkg.add("JuliaFormatter")
```

### Julia Pkg prompt

See [Start Julia Pkg prompt](@ref).

```julia
(@v1.8) pkg> add Revise
(@v1.8) pkg> add JuliaFormatter
```

## Open Dionysos

In VSCode, do `File/Open Folder.../` and select the folder `.julia/dev/Dionysos` inside your home directory.
Before doing any changes, make sure to [Switch to the master branch and update it](@ref); see [Workflow](@ref).

## Build the documentation

To build the documentation, start by activating the documentation environment and using the Dionysos version in development.
Start by `]` to enter the package environment:
```julia
(@v1.7) pkg> activate docs
  Activating project at `~/.julia/dev/Dionysos/docs`

(docs) pkg> dev .
   Resolving package versions...
  No Changes to `~/.julia/dev/Dionysos/docs/Project.toml`
  No Changes to `~/.julia/dev/Dionysos/docs/Manifest.toml`
```

Once in a while you can also update with
```julia
(docs) pkg> up
    Updating registry at `~/.julia/registries/General`
    Updating git-repo `https://github.com/JuliaRegistries/General.git`
  No Changes to `~/.julia/dev/Dionysos/docs/Project.toml`
  No Changes to `~/.julia/dev/Dionysos/docs/Manifest.toml`
```

If you plan to change the documentation, it might be a good idea to use `Revise` (see [Install the Revise.jl and JuliaFormatter.jl packages](@ref)):
```julia
julia> using Revise
```

If you don't plan to test the examples, comment out the Literate part in `docs/make.jl`:
```julila
 20 # for example in EXAMPLES_SOLVERS
 21 #     literate_actions(joinpath(EXAMPLES_SOLVERS_DIR, example), OUTPUT_DIR)
 22 # end
 23 # for example in EXAMPLES_UTILS
 24 #     literate_actions(joinpath(EXAMPLES_UTILS_DIR, example), OUTPUT_DIR)
 25 # end
 26 # literate_actions(joinpath(@__DIR__, "src", "Getting Started.jl"), OUTPUT_DIR)
```
This will speed up building the documentation quite a lot.

Now, build the documentation with:
```julia
julia> include("docs/make.jl")
```

To view it, open the file `docs/build/index.html` with your web browser.
