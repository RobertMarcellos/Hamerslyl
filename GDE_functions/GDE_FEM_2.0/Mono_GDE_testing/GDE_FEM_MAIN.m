% This is the master code for the reproducing the results presented in the
% article Application of the finite element method to General Dynamic
% Equation of Aerosols (GDE-FEM 1.0) - comparison with classical numerical
% approximations. Default is that this code can be run to produce all
% results presented in the article. 
close all
clear all
clc

% Add functions subfolder into the search path
addpath(genpath([fileparts(pwd),'\GDE_functions']))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FE marix creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decision for the test case FE matrix creation. Set corresponding cell 
% as 1 in order to create FE matrix

% 1 Coagulation (constant kernel)
% 2 GDE (condensation and coagulation)
% 3 discrete GDE
% 4 Condensation
create = [1 1 1 1];

% FE matrix dimensions are different for the different test cases. Set here
% the wanted sizes for different test cases.
sizes = [25,50,75,100,125,150,200,250,300,350,400,450,500,600,700,...
    800,900,1000]';
sizes_cond = [sizes; 1500; 2000; 2500; 3000; 3500; 4000; 4500; 5000];
sizes_discrete = [20,30,40,50,75,100,125,150,200,250,300,350,400,450,...
    500,600,700,800,900,1000]';


% Makes subfolder for the FE matrices
if ~exist([pwd,'\GDE_FE_matrices'],'dir')
       mkdir([pwd,'\GDE_FE_matrices'])
end
addpath([pwd,'\GDE_FE_matrices'])

FE_matrix_creator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time evolution calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Desision parameter for the time evolution computations. Set corresponding cell 
% as 1 in order to create FE matrix

% 1. Coagulation test case
% 2. Condensation test case
% 3. Analytical GDE comparison test case
% 4. Discrete GDE test case (Another decision parameter for this below)

time_evol_calc = [1 1 1 1]';

% Makes subfolder for the GDE time evolutions
if ~exist([pwd,'\Time_evolutions'],'dir')
       mkdir([pwd,'\Time_evolutions'])
end
addpath([pwd,'\Time_evolutions'])


% Constant coagulation time evolution
if time_evol_calc(1) == 1
    constant_coagulation_time_evolution
end

% Condensation Equation test case
if time_evol_calc(2) == 1
    condensation_time_evolution
end

% GDE test cases
if time_evol_calc(3) == 1
    GDE_time_evolution
end

% Discrete GDE test case
if time_evol_calc(4) == 1
    
    % Decision for creating the discrete GDE (computationally heavy) y/n
    create_discrete = 'y';
    
    discrete_GDE_time_evolution
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decide plots
% 1. condensation
% 2. coagulation
% 3. GDE analytical
% 4. discrete GDE
plot_case = [1 1 1 1];

% Choose the discretization for the figures
discretization = [200,100,100,30];

% Save format decision for the subplots
save_format = 'png';
% save_format = 'epsc';

% Creates directories for the figures of they do not exist
if ~exist([pwd,'\figs'],'dir')
       mkdir([pwd,'\figs'])
end
if ~exist([pwd,'\figs\subplots'],'dir')
       mkdir([pwd,'\figs\subplots'])
end

% Save location for the subplots
sl_sub = [pwd,'\figs\subplots\'];

% Makes plots according to set decision
plot_GDE_figs

% Example of the basis functions
test_function_plot



