# Forking and Cloning tRIBS
If you intend to develop, fix, or modify tRIBS source code it's important that you:
1) Start by creating a [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) of the tRIBS repository.
2) Create a development branch to add any changes in the source code.
3) Keep your main and development branch up to date with changes from the main tRIBS repository.
4) Updates the main tRIBS source code from your local fork can by creating [pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) on GitHub.

## Getting started with GitHub
Any modification or fixes of the tRIBS source code should take place in your own *fork* of the main tRIBS repository. A fork is a *mirror* of the repository and is hosted on your
personal GitHub account. Any updates, modifications, or fixes to the tRIBS source code can be implemented on this fork. These updates then can be merged to the tRIBS main repository through implementing a pull request on the GitHub website. It is important to note that you will need to 
create a GitHub account and use Git in order work with tRIBS source code in this way. Here are links for creating a [GitHub account](https://github.com) and [installing Git](https://help.github.com/en/github/getting-started-with-github/set-up-git), or alternatively you can use a GitHub graphical user interface [here](https://desktop.github.com). If you use Git through the command line
you will need to configure your account for write access using an SSH key with relevant documentation provided [here](https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh).

## Forking tRIBS

The following steps will create a fork of the tRIBS repository under
your GitHub account.

1. Sign in to your GitHub account.
2. Go to the tRIBS [home page]()
   on GitHub.
3. Click the fork button in the upper-right corner of the page.

Once completed, you will be redirected to the home page for your own
copy of tRIBS.

## Cloning your fork to your computer
You can clone the fork (that lives on the GitHub website) locally to
your computer either using the GitHub app (the GUI), or directly from
the command line using git. If you've never used git before, the app is
probably the way to go. 

### Using the GUI/app

1. Sign in to GitHub on the GUI.
2. Click on your account on the left side of the GUI.
3. This will show you your fork of tRIBS. Click on the clone option next to
   the fork which will initiate the download of tRIBS to your computer.

### Using the command line

Use the following commands from the terminal.

``` bash
   $ git clone git@github.com:your-user-name/tribs_sub2020.git
   $ cd tribs_sub2020
   $ git remote add upstream git://github.com/tribshms/tribs_sub2020.git
```
Note if you are running into the following errors,
``` bash
git@github.com: Permission denied (publickey).
fatal: Could not read from remote repository.
```
then you will likely need to generate a new ssh key and add it to your [GitHub profile](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).