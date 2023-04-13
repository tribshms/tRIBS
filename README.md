# tRIBS 2023
tRIBS source files for development 

This version of tRIBS has been updated with a [CMake](https://cmake.org/) build system and has been updated to currently run on new Mac OS (2022-2023) with C++ standard 17.

Click [here](./md/CMake.md) for instructions using CMake to compile tRIBS.

To see the latest changes incorporated into this version of tRIBS see the [CHANGELOG](./md/CHANGELOG.md) and for information on upcoming updates and tasks for tRIBS development see [TODO](./md/TODO).


<!-- This verison of tRIBS incorporates significant developement, with the latest modifications made by Josh Cederstrom over the course of his degree. This version also includes fixes from Ara Ko, Carlos Lizarraga, and Xiaoyang. It does not include updates to make files, as a current goal is to bypass this step by using CMake.

Fixes includes:
- various bug fixes that may or may have not been incorporated into the main code version
- additions to model outputs to meet the needs of my project e.g. adding snowpack sublimation and evaporation to the outputs
- changes in variable outputs for certain files e.g. replacing certain variables in dynamic file with the cumulative values
- additions to snow model like windspeed reduction below the canopy and ability to specify snow liquid water holding capacity in the input file.

To find these changes search for "CJC", "CJC2019", "CJC2020", and "CJC2021"

This version does incorporate the ability to load soil parameters in the form of grids. -->
