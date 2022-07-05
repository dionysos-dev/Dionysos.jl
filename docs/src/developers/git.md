# Git recommended workflow

Git is very flexible and this can be a bit too much at first.
This guide provides a workflow that should allow you to get things done and not lead you in any tricky situations.
Moreover, it guarantees that the whole team has the same setup.

## Add your remote

### VSCode

In VSCode, you should do `Shift+Ctrl+G` then click on the three dots at the top right of the left pane then `Remote` then `Add remote...` then enter `https://github.com/jdupont/Dionysos.jl.git` and then `jdupont`.

### Command line

Suppose your Github login is `jdupont`, add your fork (assuming you have already clicked on the "Fork" button on Github) as a remote:
```sh
$ git remote add jdupont https://github.com/jdupont/Dionysos.jl.git
```
Your remotes should be (the order of the lines is not important):
```
$ git remote -v
jdupont	https://github.com/jdupont/Dionysos.jl.git (fetch)
jdupont	https://github.com/jdupont/Dionysos.jl.git (push)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (fetch)
origin	https://github.com/dionysos-dev/Dionysos.jl.git (push)
```

## Switching to the master branch and updating it

### VSCode

Click on the branch on the lower left and enter `master`.
Now, click on the rotating arrows on the lower left at the right of `master` to update it.

### Command line

Before you start working on something new, pull any new changes made by the team to the master branch of your computer.
```sh
$ git checkout master # Switch to the master branch of your computer
$ git fetch origin master # Fetch the new commits of the master branch on Github
$ git merge --ff-only origin/master # Merge the new commits into the master branch of your computer
```

## Creating a new branch

Start by [Switching to the master branch and updating it](@ref).
Create a new branch (choose a branch name, let's suppose it is `mybranch`) and switch to it.

### VSCode

Click on `master` on the lower left and then `+ Create new branch...` then write `mybranch`.

### Command line

```sh
$ git branch mybranch # Creates a new branch `mybranch`
$ git checkout mybranch # Switch to the new branch `mybranch`
```

## Pushing changes

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

## Peer reviewing

Now go on your fork on the Github website and open a pull request.
You should receive reviews asking you to do changes.
Do these changes on your computer and push them as explained [above](pushing-changes).
Once your changes are accepted and merged, delete your branch with `git branch -d` as follows.
Start by [Updating the master branch](#updating-the-master-branch),
otherwise, `git branch -d` won't see that the branch is merged and will disallow to delete it to avoid losing your work.
```sh
$ git checkout master
$ git branch -d mybranch
```
Now go back to [Creating a new branch](#creating-a-new-branch), it's easier to use a new branch instead of using the same branch again.

## Resolving conflicts

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
