These MATLAB codes can be used to simulate the time evolution of aerosol number distribution undergoing condensation and coagulation processes. The finite element method is applied for the approximation method,

The results presented in article "Application of finite element method to General Dynamic Equation of Aerosols (GDE-FEM 1.0) - comparison with classical numerical approximations" can be replicated by running the code
GDE_FEM_MAIN

GDE_FEM_MAIN calls for the other .m -files.

In addition, the example_size_dist_evolution.m is included which can be used to test the different functions included in this package.


FUNCTIONS EXPLAINED (Inside GDE_functions folder)

------------------------------------------- FE matrix creation codes -----------------------------------------------

Coagulation_FEM_matrix_creator: Creates coagulation FEM matrices with trapezoidal rule method

Coagulation_quadrature_matrix_creator: Creates coagulation FEM matrices with 3 or 5 point Gaussian quadrature

Condensation_FEM_matrix_creator: Creates Petrov-Galerkin and Galerkin FEM matrixes for the condensation equation

------------------------------------------- Size splitting operator for discrete GDE -----------------------------------------------

Size_splitting_operator: Size splitting operator for discrete coagulation equation

------------------------------------------- Time evolution codes -----------------------------------------------

CrankNicolsonGDE: Crank-Nicolson mainly designed for the FEM time evolution (Should also work with difference method).

discrete_GDE_solver_ver2: Calculates solution for discrete General dynamic equation. Monomer size have to be set and a minimum particle size. 
			  Can be used to calculate only coagulation if source is set to zero.

------------------------------------------- Plotting codes -----------------------------------------------

fig: Code for configuring figures

PlotParticleDensityEvolution: Colorplot with logarithmic axices
