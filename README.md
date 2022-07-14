# cortex_patterning
Repository of codes related to simulating the mathematical model and performing statistical analysis for Feng et al., Science Advances, 2021 (see https://doi.org/10.1126/sciadv.abf6808).  

The collection of code here is for:

1. solving PDEs and processing output for models of cortical patterning/neural fate in brain development.
2. Processing cortical imaging data and computing statistics on this data - see the imaging_analysis folder.

Files with the .edp extension are for the freely available finite element PDE solver FreeFEM (https://freefem.org/). 

* The ffmatlib folder contains codes needed to interpret FreeFEM meshes and solution output in MATLAB.
* ffmatlib.idp is a library which contains functions to write FreeFEM output to files (it should be in the same parent folder as the .edp file you want to run).

## 1D Simulations
These codes produce Math Supplement Figure 6 from Feng et al. (2021)
* 1D_gradient_linear.edp.edp
* matlab_plotting_script.m

## 2D Simulations
These codes produce Math Supplement Figure 7 from Feng et al. (2021)
* 2D_square_linear_gradient.edp
* plot2D_square.m

## 3D Simulations
These codes produce Math Supplement Figure 10 from Feng et al. (2021)
* 3D_simulation.edp
* plot3D_slice.m 

