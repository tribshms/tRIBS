| Operating System | Build Status |
|------------------|--------------|
| Linux            | ![Linux Build](https://img.shields.io/github/actions/workflow/status/tribshms/tRIBS/compile_and_test_linux.yml) |
| macOS            | ![macOS Build](https://img.shields.io/github/actions/workflow/status/tribshms/tRIBS/compile_and_test_macos.yml)|
| Windows          | *Not Supported* |

![](https://img.shields.io/readthedocs/tribshms)

# TIN-based Real-time Integrated Basin Simulator: Version 5.2
This repository contains source code for the fully distributed hydrological model: TIN-based Real-time Integrated Basin Simulator (tRIBS). Details on running the model and its applications can be found [here](https://tribshms.readthedocs.io/en/latest/).

Licensing information can be found in [LICENSE.md](./LICENSE.md).

## Installation and Compilation
tRIBS is an object-oriented based code written in c++. In order to run tRIBS you must download and compile the source code from this GitHub repository.
Instructions for using CMake to compile tRIBS can be found [here.](./doc/md/CMake.md)

Alternatively, a Docker image of tRIBS is maintained on [Docker Hub](https://hub.docker.com/repositories/tribs). 
Further information can be found [here](doc/md/DOCKER.md).


**For additional details in working with source code, including pertinent git commands, see the following instructions [here](./doc/md/DEV_INST.md). Additional information about debugging tRIBS can be found [here](./doc/md/DEBUG.md)**
