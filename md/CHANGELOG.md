<!--- CHANGELOG.md --->
# Changelog

All notable changes to this project will be documented in this file.

<!--- 
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
--->

## [Unreleased]

## [V 3.1] - 2023-4-05

### Added
- CMake functionality (CMakeLists.txt)
- Merged fixes from different versions of tRIBS code including from Josh Cederstrom, Ara Ko, Carlos Lizarraga, and Xiaoyang Tang
- added #include "tTimer.h" to tTimer.cpp
- markdown (md) sub directory to display markdown files on github


### Fixed
- Compiler errors related to assert statements with null pointers
- Compiler error for tPtrList.h (L 873) newlist.insertAtBack to newlist->insertAtBack 
- If optres == 0, would define tReservoir/tResData as dangling pointer, fixed by making tReservoir public and setting as null if optres == 0.
- Commented out #include t*(parallel code).cpp in parallel header files because it led to redefinition

### Changed

### Removed
- removed register calls (no longer supported at c++ 17 or earlier)
- removed old make files

