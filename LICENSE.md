# TIN-based Real-time Integrated Basin Simulator (tRIBS) Distributed Hydrological Model
## Version 5.2.0

tRIBS Distributed Hydrologic Model, Version 5.2.0, represents a culmination of efforts and significant improvements from the intial release[^1] and later versions. This includes increases in computational efficiency through parallelization[^2], ingestion and use of dynamic land use grids derived from remote sensing products[^3], and incorporation of additional physical processes, ranging from a single-layer snowpack[^4], to level-pool reservoir routing[^5] and channel transmission losses[^6]. The latest updates (Version 5.0 and onwards) entail a CMake build system and major code improvements, including a refactored snow module, fixed memory leaks, and updates to C++ 17 standards. 

## MIT License

Copyright © 2024 tRIBS Developers

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## References:
[^1]: 1. Ivanov, V. Y., Vivoni, E. R., Bras, R. L. & Entekhabi, D. Catchment hydrologic response with a fully distributed triangulated irregular network model. Water Resources Research 40, (2004).
[^2]: Vivoni, E. R. et al. Real-world hydrologic assessment of a fully-distributed hydrological model in a parallel computing environment. Journal of Hydrology 409, 483–496 (2011).
[^3]: Vivoni, E. R. Diagnosing Seasonal Vegetation Impacts on Evapotranspiration and Its Partitioning at the Catchment Scale during SMEX04–NAME. Journal of Hydrometeorology 13, 1631–1638 (2012).
[^4]: Rinehart, A. J., Vivoni, E. R. & Brooks, P. D. Effects of vegetation, albedo, and solar radiation sheltering on the distribution of snow in the Valles Caldera, New Mexico. Ecohydrology 1, 253–270 (2008).
[^5]: Cázares-Rodríguez, J. E., Vivoni, E. R. & Mascaro, G. Comparison of Two Watershed Models for Addressing Stakeholder Flood Mitigation Strategies: Case Study of Hurricane Alex in Monterrey, México. Journal of Hydrologic Engineering 22, 05017018 (2017).
[^6]: Schreiner-McGraw, A. P. & Vivoni, E. R. On the Sensitivity of Hillslope Runoff and Channel Transmission Losses in Arid Piedmont Slopes. Water Resources Research 54, 4498–4518 (2018).
