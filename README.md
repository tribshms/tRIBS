# TIN-based Real-time Integrated Basin Simulator
 The TIN-based Real-time Integrated Basin Simulator (tRIBS, version 4.0) is a fully distributed physical based hydrological model. Details about the running the model and its applications can be found [here](https://tribshms.readthedocs.io/en/latest/).
## Installation and Compilation
tRIBS is an object-oriented based code written in c++. In order to run tRIBS you must download and compile the source code from this GitHub repository. tRIBS V4 has been updated with a CMake build system and has been updated to currently run on new (2022-2023) Mac and Linux OS with C++ standard 17.
Instructions for using CMake can be found [here.](./md/CMake.md)
## Working with tRIBS source code
If you intended to develop, fix, or modify tRIBS source code it's important that you:
1) Start by creating a [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) of the tRIBS repository.
2) Create a development branch to add any changes in the source code.
3) Keep your main and development branch up to date with changes from the main tRIBS repository.
4) Updates the main tRIBS source code from your local fork can by creating [pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) on GitHub.

For additional details in working with source code, including pertinent git commands, see the following instructions [here](./md/DEV_INST.md)
