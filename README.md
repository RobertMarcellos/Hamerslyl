# GDE-FEM_2.0

<a href="https://doi.org/10.5281/zenodo.8092361"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4769125.svg" alt="DOI"></a>

MATLAB codes for simulating the time evolution of aerosol number distribution undergoing condensation, coagulation, removal and source processes with Finite Element Method (FEM) and Petrov-Galerkin Finite Element Method (PGFEM).
The aerosol particles can contain multiple compounds. Temporal evolutions are calculated with the Crank-Nicolson method. Code package includes support for FEM approximation for single and multicomponent GDE. 
Codes are individually commented and this documents gives an overview of the functions.


```sh
example_size_dist_evolution.m
```
The example_size_dist_evolution.m uses the functions provided in this package.


The results presented in article "Teemu Salminen, Kari E.J. Lehtinen, Jari P. Kaipio, Vincent Russell, Aku Seppänen, Application of finite element method to General Dynamic Equation of Aerosols (GDE-FEM 1.0) - comparison with classical numerical approximations, 
Journal of Aerosol Science, https://doi.org/10.1016/j.jaerosci.2021.105902" can be replicated by running  code in folder "Mono_GDE_testing":

```sh
GDE_FEM_MAIN.m
```

The results presented in article "Application of finite element method to multicomponent General Dynamic Equation of Aerosols" can be replicated by running the code by running code in folder "MC_GDE_testing":
```sh
MC_GDE_MAIN.m
```

## Provided functions
Functions for the aerosol number distribution are located in the subfolder GDE_functions. Codes for the multicomponent GDE are in subfolder `multi_component`, and codes for the single component are in subfolder  `single_component `.
Short explanation for the codes are given here.

### FE matrix creation codes

```sh
Single component
```

- `Coagulation_FEM_matrix_creator` : Creates coagulation FEM matrices with trapezoidal rule method.

- `Coagulation_quadrature_matrix_creator` : Creates coagulation FEM matrices with 3 or 5 point Gaussian quadrature.

- `Coagulation_quadrature_matrix_creator2` : !!Improved version of previous code!! Creates coagulation FEM matrices with 3 or 5 point Gaussian quadrature.

- `Condensation_FEM_matrix_creator` : Creates Petrov-Galerkin and Galerkin FEM matrixes for the condensation equation.

```sh
Multicomponent
```

- `Coagulation_MC_quadrature_matrix_creator` : Creates FEM and PGFEM coagulation matrices for the approximation of MCGDE

- `multicomponent_GDE_FE_matrix_creator` : Creates FEM and PGFEM matrices for the condensation and removal process

- `multicomponent_GDE_sectional_matrices_volume` : Create matrices for temporal evolution calculated with sectional method

### Size splitting operator for sectional method approximation of GDE

- `Size_splitting_operator` : Size splitting operator for discrete coagulation equation.

### Time evolution codes 

```sh
Single component
```

- `CrankNicolsonGDE` : Crank-Nicolson mainly designed for the FEM time evolution (Should also work with difference method).

- `CrankNicolsonGDE_new` : Crank-Nicolson mainly designed for the FEM time evolution (Should also work with difference method). Implementation for the boundary conditions have been modified.  

- `discrete_GDE_solver_ver2` : Calculates solution for discrete General dynamic equation. Monomer size have to be set and a minimum particle size. 
			  	Can be used to calculate only coagulation if source is set to zero.

- `discrete_GDE_solver_final` : Calculates solution for discrete General dynamic equation. Monomer size have to be set and a minimum particle size. 
			  	Can be used to calculate only coagulation if source is set to zero. Compared to `discrete_GDE_solver_ver2` boundary conditions are applied better in this code. 

- `TimeEvolutionGDE` : Calculates time evolution for the aerosol number distribution with chosen time integration method. Possible choises are: explicit and implicit Euler method and the Crank-Nicolson method.

- `Error_estimator` : Calculates the relative error for the approximation method compared to accurate solution.

```sh
Multicomponent
```

- `multicomponent_GDE_CrankNicolson` : Calculates FEM/PGFEM approximations for temporal evolution of volume concentrations with C-N method. Condensation, coagulation, removal and source processes can be included.

- `MC_GDE_sectional_method_volume_conc` : Calculates sectional method approximations for temporal evolution of volume concentrations with C-N method. Condensation, coagulation, removal and source processes can be included.

- `discrete_multivolume_GDE` : Calculates temporal evolution of volume concentrations with multivolume discrete GDE presented in Application of finite element method to multicomponent General Dynamic Equation of Aerosols, 2023, Salminen et. al



### Plotting codes

- `fig` : Code for configuring figures.

- `PlotParticleDensityEvolution` : Colorplot with logarithmic axices.


This code package is under MIT-license, and correct citation should be used when codes from this package are used.
