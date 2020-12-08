# cortex_patterning
Repository of codes to: 

1. solve PDEs and process output for models of cortical patterning/neural fate in brain development.
2. Process cortical imaging data and compute statistics on this data

Files with the .edp extension are for the freely available finite element PDE solver FreeFEM (https://freefem.org/). 

* The ffmatlib folder contains codes needed to interpret FreeFEM meshes and solution output in MATLAB.
* ffmatlib.idp is a library which contains functions to write FreeFEM output to files (it should be in the same parent folder as the .edp file you want to run).

## 1D Simulations

## 2D Simulations
These codes produce Math Supplement Figure 7
* 2D_square_linear_gradient.edp
* plot2D_square

## 3D Simulations
These codes produce Math Supplement Figure 10
* 3D_simulation.edp
* plot3D_slice.m 

