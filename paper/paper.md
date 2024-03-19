---
title: 'tRIBS v5.2: A multi-resolution, parallel platform for tributary hydrology in forest applications'
tags:
- C++
- hydrology
- distributed hydrological models
- numerical modeling

authors:
- name: Raming, L. W.
orcid: 0000-0003-2204-4621
affiliation: "1, 2" # (Multiple affiliations must be quoted)
corresponding: true # (This is how to denote the corresponding author)

- name: Vivoni, E. R.
orcid: 0000-0002-2659-9459
affiliation:  "1, 2"

- name: G.,Mascaro
orcid:0000-0003-4516-1206
affiliation:  "1, 2"

- name: Cederstrom, C.J.
affiliation: 1

-name: Ko, A.
affiliation: 1

-name: Schreiner-McGraw, A.P
affiliation: 1
orcid: 0000-0003-3424-9202

-name: Cazares-Rodriguez, J.E.
affiliation: 1

-name:  Lizarraga-Celaya, C.
affiliation: 1
orcid: 0000-0002-0893-4268


affiliations:
index: 1
- name: School of Sustainable Engineering and the Built Environment, Arizona State University, Tempe, AZ, USA, 85287.
index: 2
- name: Center for Hydrologic Innovations, Arizona State University, Tempe, AZ, USA, 85287.

date: 19 March 2024
bibliography: paper.bib
---

# Summary
Distributed hydrologic models provide earth scientists and engineers with tools to test and explore hypotheses related to the movement and storage of water within a landscape `[@fatichi2016; @grayson2002; @keller2023]`. The Triangulated Irregular Network (TIN)-based Real-Time Integrated Basin Simulator (tRIBS; `@ivanov2004a`, `@ivanov2004b`), an example of such a process-based distributed model, has been used to address a wide range of problems from hillslope scale processes in ecohydrology (e.g. `@mahmood2011`) to flood management of large watersheds (e.g. Cázares-Rodríguez et al., 2017). Yet, in spite of the extensive use and application of tRIBS to current topics in hydrology, engineering, and the earth sciences, the code has been essentially maintained as a proprietary software. Here, we document the release of tRIBS v5.2, an updated open source code base and its application for forested watersheds that serve as tributaries to larger river systems. This release includes improvements in hydrologic processes (i.e. channel transmission losses and reservoir routing; `@schreiner-mcgraw2018`; `@cazares-rodriguez2017`) as well as an updated documentation, improved infrastructure for sustainable code development and employment, and improved computational efficiency. These additions provide a robust and sustainable code base that will enhance access and applications of the model.
# Statement of Needs
## Model Description
tRIBS is written in C++ and uses an object oriented design founded on a hydrologically conditioned triangulated irregular network (TIN) mesh (see Figure 1; `@tucker2001`; Vivoni et al., 2004). Building on the work of `@garrote1995`, tRIBS is a continuous hydrologic model simulating the coupled dynamics between the vadose and saturated zones `[@vivoni2007]`. Accounting for these key hydrologic processes while using computationally efficient methods \ref{Figure 1.}, tRIBS actively tracks both the evolution of wetting fronts and moisture losses, allowing for continuous simulation throughout wet and dry periods `[@ivanov2004a]`. With the addition of a single-layer snowpack module `[@rinehart2008]`, tRIBS can also be applied in cold and mountainous forest environments \ref{Figure 2.}.
Furthermore, the unstructured TIN mesh in tRIBS provides a multiresolution approach to distributed hydrologic modeling (\ref{Figure 1.} & \ref{Figure 3.}). As a consequence, tRIBS allows for detailed control in resolving hydrologic dynamics across multiple scales `[@vivoni2004]`, maximizing model fidelity to physical processes, while minimizing computational expenses. This multi-scaling behavior when paired with parallelization `[@vivoni2011]` allows for hyper-resolution modeling `[@wood2011]` of hydrologic dynamics at unprecedented scales, from simulations rendered at a point (`@vivoni2010`, \ref{Figure 2.}) to 21,000 km2 watersheds simulated for a period of 10 years at a nominal cell resolution of ~88 m `[@ko2019]`.
![Figure 1. Conceptual overview of tRIBS end-to-end workflow highlighting key processes. Asterisks indicate new features or processes available in tRIBS v5.2. Soil and vegetation parameters may be provided in a raster with continuous values or in a classification table.[]{label="Figure 1."}](figures/Fig_1.png)
## Updates and Modifications
Building on core tRIBS functionality, as described by `@ivanov2004a`, `@rinehart2008`, and `@vivoni2011`, tRIBS v5.2 provides two new process additions: (1) reservoir routing using the level-pool method `[@cazares-rodriguez2017]`, and (2) channel transmission losses `[@schreiner-mcgraw2018]`. In addition, the code base has been restructured with mechanisms for improved maintainability, robustness, performance, and integration. This includes updates for code compatibility with newer compilers (Clang and GCC), the introduction of a CMake build system providing flexibility for compiling serial and parallel versions, and the modernization of the model version control system and documentation. Additionally, we refactored the snow module `[@rinehart2008]`, resulting in a reduction of redundant code and enhanced code organization. Memory leaks associated with parallel operations were fixed, allowing for increased scalability. Finally, we included Docker images for both tRIBS v5.2 and the auxiliary program MeshBuilder. The Docker image for MeshBuilder facilitates an end-to-end workflow that utilizes METIS `[@Karypis1998]`, enabling rapid and easy partitioning of a watershed domain for parallel simulations `[@vivoni2011]`.
These and other features of tRIBS v5.2 can be explored using two newly updated benchmark scenarios. This first benchmark is a point-scale simulation of the Happy Jack SNOTEL site in northern Arizona, USA (\ref{Figure 2.}). The second is a basin-scale simulation of the Big Spring watershed located in the headwaters of Sycamore Creek in northern Arizona (\ref{Figure 3.}). Both benchmarks are hosted on Zenodo, see \ref{Figure 2.} and \ref{Figure 3.}for more details.
![Figure 2. A point-scale (i.e. a single Voronoi cell) tRIBS simulation of snow water equivalent (SWE) at the Happy Jack SNOTEL site in northern Arizona, USA. Top panel shows the time series of observed (black) and simulated SWE (blue). Bottom panel compares the observed and simulated peak SWE from 2002 to 2017. Dashed black line is a one-to-one relation. The color bar indicates the time difference in the occurrence of the peak SWE for each water year. Zenodo repository for this simulation with additional details can be found at: NEED TO UPDATE []{label="Figure 2."}](figures/Fig_2.png)
![Figure 3. An example of a basin-scale tRIBS simulation showing a spatial map of mean hourly evapotranspiration rates averaged over the course of a 4-year simulation period. Big Spring basin is a tributary to Sycamore Creek in northern Arizona, USA. Zenodo repository for this simulation with additional details can be found at: NEED TO UPDATE[]{label="Figure 3."}](figures/Fig_3.png)
## Conclusion
Embracing the FAIR principles (Findability, Accessibility, Interoperability, and Reusability; `[@wilkinson2016`) and recognizing the importance of free and open source software in hydrology `[@Kabo-bah2012]`, here we document the release of tRIBS v5.2. This version represents years of cumulative efforts with major code improvements related to maintainability, robustness, performance, and integration as well as new process based functionality. The model applications included as benchmark cases serve as examples of forest applications in tributary watersheds to larger river systems. We anticipate that tRIBS v5.2 will be a valuable asset in addressing a wide range of problems for the broader hydrology community.

# Acknowledgements
We thank the tRIBS developers at the Massachusetts Institute of Technology, New Mexico Institute of Mining and Technology, Los Alamos National Laboratory, and Arizona State University, including Valeriy Y. Ivanov, Scott M. Rybarczyk, Greg E. Tucker, Sue Mniszewski, Pat Fasel, and Alex J. Rinehart. We also thank Rafael L. Bras, Dara Entekhabi, and Everett P. Springer, for guidance on model development. Lastly, we are grateful for the support of Elvy Barton and Bruce Hallin for encouraging further development and application of tRIBS to new and pressing problems. Over the years, tRIBS model development has been funded by: Army Research Office, National Science Foundation, National Oceanic and Atmospheric Administration, National Aeronautics and Space Administration, Los Alamos National Laboratory, Salt River Project, and Arizona State University. The most current contribution documented here was facilitated by the Arizona Water Innovation Initiative.




