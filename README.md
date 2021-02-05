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

## Git recommended workflow

Git is very flexible and this can be a bit too much at first.
This guide provides a workflow that should allow you to get things done and not lead you in any tricky situations.

### Set up

First, clone Dionysos:
```sh
$ git clone https://github.com/dionysos-dev/Dionysos.jl.git
$ cd Dionysos.jl
```
Suppose your Github login is `jdupont`, add your fork (assuming you have already clicked on the "Fork" button on Github) as a remote
```sh
$ git remote add jdupont https://github.com/jdupont/Dionysos.jl.git
```
Your remotes should be (the order of the lines is not important)
```
$ git remote -v
jdupont	https://github.com/jdupont/Dionysos.jl.git (fetch)
jdupont	https://github.com/jdupont/Dionysos.jl.git (push)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (fetch)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (push)
```

### Updating the master branch

Before you start working on something new, pull any new changes made by the team to the master branch of your computer.
```sh
$ git checkout master # Switch to the master branch of your computer
$ git fetch origin master # Fetch the new commits of the master branch on Github
$ git merge --ff-only origin/master # Merge the new commits into the master branch of your computer
```

### Creating a new branch

Create a new branch (choose a branch name, let's suppose it is `mybranch`) and switch to it with:
```sh
$ git branch mybranch # Creates a new branch `mybranch`
$ git checkout mybranch # Switch to the new branch `mybranch`
```

### Pushing changes

Before doing any changes, check that your are on the right branch.
```sh
$ git checkout mybranch # Switch to the branch `mybranch`
```
Now, do you changes...

Once you are done, check what you want to commit with `git status`.
For instance, if you want to commit every change in the `src` folder in addition to the changes in `test/sometest.jl` but not the rest, do `git add src` and `git add test/sometests.jl`.
If you want to commit everything, do `git add .`.
Check what you have added with `git status`.
If you want to remove what you have added, do `git reset`.
Once `git status` shows the correct output, do
```sh
$ git commit -m "Short summary of what you have done"
$ git push jdupont mybranch
```

### Peer reviewing

Now go on your fork on the Github website and open a pull request.
You should receive reviews asking you to do changes.
Do these changes on your computer and push them as explained [above](pushing-changes).
Once your changes are accepted and merged, start by [Updating the master branch](#updating-the-master-branch),
otherwise, `git branch -d` won't see that the branch is merged and will disallow to delete it to avoid losing your work.
```sh
$ git checkout master
$ git branch -d mybranch
```
Now go back to [Creating a new branch](#creating-a-new-branch), it's easier to use a new branch instead of using the same branch again.

### Resolving conflicts

By the time you have created new branch, other developers may have made changes to the master branch at the same lines of the same files as the changes in some of your branch.
In that case, Github won't allow your pull request to be merged as it does not know whether it should take the changes of your branch or the changes of master.
To fix the conflicts, start [Updating the master branch](#updating-the-master-branch) and then do
```sh
$ git checkout mybranch
$ git rebase master
```
The rebase will stop at every conflicting commit and you will have to, edit the conflicting files and edit parts of the file like
```
<<<<<<< HEAD
changes on master
=======
changes on your branch
>>>>>>> mybranch
```
and choose the changes that should be kept. Then you should `git add` the conflicting files, run `git rebase --continue`.
See [here](https://docs.github.com/en/github/using-git/resolving-merge-conflicts-after-a-git-rebase) for more information on this.
In case you have may conflicting commits, the rebase can be tedious, you will even have to resolve conflicts for some changes in some earlier
commit even if you reverted these changes in a later commit. To avoid this issue, it is recommended to first squash your commits into a single one.
To do that, use
```sh
$ git rebase -i master
```
and replace `pick` by `s` for all commits except the top one and then follow the instructions, see [here](https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History)
for more details.
