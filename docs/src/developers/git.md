# Git recommended workflow

Git is very flexible and this can be a bit too much at first.
This guide provides a workflow that should allow you to get things done and not lead you in any tricky situations.
Moreover, it guarantees that the whole team has the same setup.
This guide assumes you have already followed and completed your [Set up](@ref).
This guide assumes that your Github login is `jdupont`. Replace it by the appropriate login since that is not the case.

You can be in 2 situations:
1) You don't have write access to `origin`. If you don't know what this means, it means you are in that situation.
   Once you have enough mastery of Git, we might give you write access but it's best to start in this situation.
2) You have write access to `origin`. This means you can push directly to `master` which you should never do. We only give
   you write access when we trust that you won't mess up.

We have a `master` branch that contains the latest version of all **merged** changes.
There is three `master` branches:
* `origin/master`: that is the branch at https://github.com/dionysos-dev/Dionysos.jl, it is always the most up to date.
* Your local `master` branch: that is the state of the branch on your computer. It may be a few commits behind `origin/master` as they do not synchronize automatically. You can update it by following [Switch to the master branch and update it](@ref).
* `jdupont/master`: that is the branch of your fork https://github.com/jdupont/Dionysos.jl, it may be many commits behind your local `master` and even more commits behind `origin/master` but we don't care much because we won't use it so you don't have to update it. Once you have write access to `origin`, you care even less if that was possible.

## Workflow

The workflow is as follows. Your contributions should be grouped into small chunks that bring Dionysos from a working state (which is the current version of `origin/master`) to a new working state containing your improvements.
In order to do that, you do the following for every chunk.
1) [Switch to the master branch and update it](@ref);
2) then you [Create a new branch](@ref) (let's call it `mybranch` but it should be a new name for every new small chunk);
3) then you make your changes on your computer...
4) then you [Format your code](@ref) with [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl), so that the whole code has the same format; 
5) then you [Commit your changes](@ref), this will update your local version of `mybranch`;
6) then you [Push your changes](@ref), this will update the version of `mybranch` in your fork (resp. `origin`) if you don't have write access) (resp. you have write access)) ; 
7) then you [Create a pull request on Github](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request), note that if the code was not properly formatted, the pull request will raise an error on Github;
8) then you should receive reviews asking you to do changes. Do these changes on your computer and push them. To do that, [Switch branches](@ref) to `mybranch` then follow steps 4) and then 5) again;
9) you might need to [Resolve conflicts](@ref), especially if you did not follow step 1);
10) Once your branch has been merged, [Delete your branch](@ref) and go back to step 1) for a new chunk.

If you want to make a change independent from the change you are currently making go back to step 1) use another branch name, say `myotherbranch`.
Once you want to go back to the changes you were doing previously, [Switch branches](@ref) back to `mybranch`.
Changing branches requires you that you have committed your changes so you should at least finish step 4) in order to do that.
It's best to also have completed step 5) so that your changes are backed up in the cloud in case something happens to your computer.

## Fork and add your remote

These are the steps that should be done only the first time when you set up.

At first, you don't have write access, so create a fork https://github.com/jdupont/Dionysos.jl by going to https://github.com/dionysos-dev/Dionysos.jl and click on the "Fork" button on the top right.
You should then add this as a remote as detailed below.

### VSCode

Switch to Source Control by pressing `Ctrl+Shift+G` then on the three horizontal dots on the top right of the left pane then `Remote` then `Add remote...` then enter `https://github.com/jdupont/Dionysos.jl.git` (replace `jdupont` by your Github login!) and then `jdupont`.

### Git bash

First [Start Git bash](@ref).

Then, add the remote as follows:
```sh
$ git remote add jdupont https://github.com/jdupont/Dionysos.jl.git
```
Your remotes should be (the order of the lines is not important):
```sh
$ git remote -v
jdupont	https://github.com/jdupont/Dionysos.jl.git (fetch)
jdupont	https://github.com/jdupont/Dionysos.jl.git (push)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (fetch)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (push)
```

## Switch to the master branch and update it

This should be before any new change! See [Workflow](@ref).

If the procedure below fails or you get a message about a need to create a merge, it means you have commit changes to your local master branch, you did not follow the [Workflow](@ref).
Contact us to get help, we will be mad at you for not following the [Workflow](@ref) but we will still help you.

### VSCode

Click on the branch on the lower left and enter `master`.
Now, click on the rotating arrows on the lower left at the right of `master` to update it.

### Git bash

First [Start Git bash](@ref).

Before you start working on something new, pull any new changes made by the team to the master branch of your computer.
```sh
$ git checkout master # Switch to the master branch of your computer
$ git fetch origin master # Fetch the new commits of the master branch on Github
$ git merge --ff-only origin/master # Merge the new commits into the master branch of your computer
```

## Create a new branch

Start by [Switch to the master branch and update it](@ref).
Create a new branch (choose a branch name, let's suppose it is `mybranch`) and switch to it.

### VSCode

Click on `master` on the lower left and then `+ Create new branch...` then write `mybranch`.

### Git bash

First [Start Git bash](@ref).

```sh
$ git branch mybranch # Creates a new branch `mybranch`
$ git checkout mybranch # Switch to the new branch `mybranch`
```

## Commit your changes

Before doing any changes, make sure you [Switch branches](@ref) to the right branch (which should be `mybranch`).
Once you have made changes, they are saved on your disc but the Git history has not been modified yet; neither your local branches nor the remote one on `origin` or your fork!
You should first stage changes (this means selecting changes you want to commit and hence to be applied on the chunk of changes you suggest to make to Dionysos in the Pull Request) and then create the commit with a message.
Usually, you want to stage every file; both the *modified* ones and the *untracked* ones because the "untracked" ones you don't want to add should be listed in the `.gitignore` file and hence should not show up.
The only exception is the file `docs/Project.toml` since you may have changed it by adding `Dionysos` in [Build the documentation](@ref) but you don't want to push that, you want to keep these changes to your computer only.

### VSCode

Switch to Source Control by pressing `Ctrl+Shift+G`.

You will see the modified files with a `M` and the untracked files with a `U`.
Click on the "+" to stage all the changes of a file.
Once all files have been staged, Write a message in the field above the "Commit" blue button.
Then press on the "Commit" blue button.

### Git bash

First [Start Git bash](@ref).

```sh
$ git checkout mybranch # Make sure you are not on the `master` branch!
$ git status # Shows modified and untracked files
$ git add foobar.jl # Stage file `foobar.jl`, replace it by the files you want to stage
$ git commit -m "Commit message" # Replace "Commit message" by a very short message about your changes
```

To unstage every file, do
```sh
$ git reset
```

Alternatively, the following adds all *modified* and *untracked* files.
```sh
$ git add .
```

Another option is the following which adds all *modified* files and commits them directly.
```sh
$ git commit -am "Commit message"
```

## Push your changes

### VSCode

Switch to Source Control by pressing `Ctrl+Shift+G` then on the three horizontal dots on the top right of the left pane then on "push".
If it is the first time you push this branch, it will ask "The branch `mybranch` has no remote branch. Would you like to publish this branch ?"; answer with "Ok".
Then it will ask which remote to push to, select `jdupont` if you don't have write access or `origin` if you have write access.

### Git bash

First [Start Git bash](@ref).

If you don't have write access, push to your fork:
```sh
$ git push jdupont mybranch
```
Otherwise, push to `origin`:
```sh
$ git push origin mybranch
```
If this fails, it means you don't have write access or that you chose a name of branch that already exists.

## Delete your branch

First, [Switch to the master branch and update it](@ref).
Otherwise, Git won't see that the branch is merged and will disallow you to delete it to avoid losing your work.

### VSCode

Switch to Source Control by pressing `Ctrl+Shift+G` then on the three horizontal dots on the top right of the left pane then on "Branch" then "Delete Branch..." then write or select `mybranch`.

### Git bash

First [Start Git bash](@ref).

```sh
$ git checkout master
$ git branch -d mybranch
```

## Resolve conflicts

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

## Switch branches

To switch to the branch `mybranch`, do the following.

### VSCode

On the bottom left, you should see a sort of "Y" symbol with empty circles at the three leaves.
On the right, you see the current branch. Click on it and then write or select `mybranch`.

### Git bash

First [Start Git bash](@ref).

```sh
$ git checkout mybranch
```

## Format your code

To format your code, run the following in your Julia REPL. Make sure you have added [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) before. 

```julia
julia> using JuliaFormatter; format(".")
```