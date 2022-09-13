# Set up

Start by installing Julia, Visual Studio Code (VSCode) and the Julia extension of VSCode as detailed [here](https://code.visualstudio.com/docs/languages/julia#_getting-started).

Launch the Julia REPL from VSCode by typing `Shift+Ctrl+P` and then `Julia: Start REPL`.
From the REPL, write `using Pkg; Pkg.develop(url="https://github.com/dionysos-dev/Dionysos.jl.git")` or `] dev https://github.com/dionysos-dev/Dionysos.jl.git`.
This will clone Dionysos with git into the folder `.julia/dev/Dionysos`.

In VSCode, do `File/Open Folder.../` and select the folder `.julia/dev/Dionysos` inside your home directory.

From the REPL, install Revise. Start by `]`:
```julia
(@v1.7) pkg> add Revise
```
We do it in the global environment so that it is available from all the environments.

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

If you plan to change the documentation, it might be a good idea to use `Revise`:
```julia
julia> using Revise
```

If you don't plan to test the examples, comment out the Literate part in `docs/make.jl`:
```julila
 11 #for example in EXAMPLES
 12 #    example_filepath = joinpath(EXAMPLES_DIR, example)
 13 #    Literate.markdown(example_filepath, OUTPUT_DIR)
 14 #    Literate.notebook(example_filepath, OUTPUT_DIR)
 15 #    Literate.script(example_filepath, OUTPUT_DIR)
 16 #end
```
This will speed up building the documentation quite a lot.

Now, build the documentation with:
```julia
julia> include("docs/make.jl")
```

To view it, open the file `docs/build/index.html` with your web browser.
