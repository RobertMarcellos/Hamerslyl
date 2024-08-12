% Master code for running results for paper Application of finite element 
% method for the multicomponent General Dynamic Equation of Aerosols,
% JOURNAL, 2023, (c) Teemu Salminen

% INITIALIZATIONS
close all
clear all
clc

% Add GDE_functions and its subfolders into search path
addpath(genpath([fileparts(pwd),'\GDE_functions']))

% Create folder for the time evolutions if it does not exist
if ~exist([pwd,'\Time_evolutions'],'dir')
        mkdir([pwd,'\Time_evolutions'])
end
addpath([pwd,'\Time_evolutions'])

if ~exist([pwd,'\figs'],'dir')
        mkdir([pwd,'\figs'])
end

% Various discretizations for the test cases. Time evolutions for the
% volume concentrations are calculated with each discretization and time
% step

discretizations = [25,50,75,100,125,150,175,200,250,300,350,400,450,500,600,700,800,900,1000,...
    1100,1200,1300,1400,1500]';


time_steps = [0.5 , 2];
time_steps_disc = [0.01, 0.04];
%%

% Multicomponent condensation test case (3 components, 2 are condensing and
% 1 is evaporating)
MC_cond_test_case

% Multicomponent condensation test case with the removal term (2 condensing
% compounds, and removal term)
MC_cond_removal_test_case

% Decision parameters for creation of the multivolume discrete GDE
% formulations. Should be done only once for set of parameters.
create_disc_coag = 'y';
create_disc_GDE = 'y';

% Discrete coagulation
if create_disc_coag == 'y'
    % NOTE: have to be run once in order to acquire accurate presentations
    % for the aerosol dynamics. 
    discrete_MC_GDE_coagulation
end
MC_discrete_GDE_coagulation_test_case

% Multivolume Discrete GDE test case (two component, Brownian coagulation kernel)
if create_disc_GDE == 'y'
    % NOTE: have to be run once in order to acquire accurate presentations
    % for the aerosol dynamics.
    discrete_MC_GDE_test_case_discrete

end
MC_discrete_GDE_test_case

%% PLOTTING CODES FOR THE TEST CASES
% save_format = 'epsc';
save_format = 'png';
% Plotting code for the test case 1.
MC_condensation_plotting
% Plotting code for the test case 2.
MC_condensation_removal_plotting
% Plotting code for the MC discrete coagulation (case 3)
MC_discrete_coagulation_plotting
% Plotting code for the MC discrete GDE (case 4)
MC_discrete_GDE_plotting
