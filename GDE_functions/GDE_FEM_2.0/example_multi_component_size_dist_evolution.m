% This code calculates an example temporal evolution with chosen parameter
% for multicomponent General Dynamic Equation of Aerosols (MCGDE), which is
% described in detail in Katoshevski, D., & Seinfeld, J. H. (1997).
% Analytical solution of the multicomponent aerosol general dynamic
% equation—without coagulation. Aerosol Science and Technology, 27(4),
% 541–549. https://doi.org/10.1080/02786829708965493
%
% Chosable parameters are:
% - Decision parameter for discrete GDE calculation
% - Plot initial distribution
% - Number of elements and bins for FEM and sectional method respectively
% - Upwinding Petrov-Galerkin parameter
% - Size interval for aerosol distribution
% - Initial distribution
% - Growth function
%
% Teemu Salminen
% Department of Technical Physics
% Unversity of Eastern Finland
% March 2023
%
clear all
close all
clc

% Add functions subfolder into the search path
addpath(genpath([pwd,'\GDE_functions']))

% CALCULATE MULTIVOLUME DISCRETE GDE
create_disc = 'y'; % 'n'/'y' (no/yes)

% PLOT INITIAL SIZE DISTRIBUTION
plot_init = 'n'; % 'n'/'y' (no/yes)

% NUMBER OF ELEMENTS AND BINS
bins = 100;
elements = 100;

% SIZE INTERVAL AS PARTICLE VOLUME
Vmin = -9;
Vmax = -6; % The size interval is chosen to be narrow because of discrete GDE

% PETROV-GALERKIN UPWINDING FACTOR
eps = 0.3;

% Discretization for the discrete GDE
vv = 10^Vmin:10^Vmin:10^Vmax;
vv = vv(:);

% Changing the value of the number distribution to value of number
% concentration (for monomers)
vN0 = vv(1);
n01 = 5e10;
N0 = vN0*n01;

% CREATING INTIAL DISTRIBUTION
sigma1 = 2e-8;
mu1 = 50*10^Vmin;
N1 = 5e6;
f1 = N1*exp(-0.5*((vv-mu1)./sigma1).^2);
n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2);
n0 = n_initial(vv);

% TIME DISCRETIZATION
tmax = 96;
delta_t = 0.5;
tt = [0:delta_t:tmax]';

if plot_init == 'y'
    figure
    semilogx(vv,n0)
    drawnow
end

% COAGULATION KERNEL FUNCTION

% Brownian coagulation kernel
% Parameters for the coagulation kernel function
T = 280; % Temperature K
rho = 1000; % kg/m^3 Density of particles
rho = rho*10^-(6*3); % kg/mum^3
k_b = 1.3806e-23; % m^2*kg*s^-2 K ^-1 Boltzmans constant
k_b = k_b*10^(2*6); % mum^2...
b = (3/(4*pi))^(1/6)*sqrt(6*k_b*T/rho);
% b = b*10^(-4*1); % unit convertion to cm/s

% Coagulation kernel function (Can be changed)
fun = @(v,w) 3600*10^(-9)*b*(1./w+1./v).^(1/2).*(v.^(1/3)+w.^(1/3)).^2;

%% Discrete GDE time evolution where smallest particle consist 2 monomers

% NOTE: Initial distributions can be set differently and condensation does
% not have to be related to number of monomers!

% Fractions for the monomers for condensation
a = 0.85;
bbb = 0.15;

% Fractions for the initial size distribution
aa = 0.5;
bb = 0.5;

if create_disc == 'y'
    [N,PHI,v_middle,v_edges,tt] ...
        = discrete_multivolume_GDE(10^Vmin,[10^Vmin,10^Vmax],@(vv) n_initial(vv),'distribution',[{aa},{bb}],...
        [a*vN0*n01,bbb*vN0*n01],fun,delta_t,tmax,'GDE');
end
%% FEM matrix creations for the problem
% Single component GDE
if create_disc == 'y'
    g = logspace(log10(2*10^(Vmin)),Vmax,elements)';
else
    g = logspace(Vmin,Vmax,elements)';
end

% SET INITIAL DISTRIBUTION FOR THE FEM and sectional (can be changed)
initial_dist = [{aa.*g.^(1).*n_initial(g)};{bb.*g.^(1).*n_initial(g)}];
initial_dist_sec = [{@(g) aa.*g.^(1).*n_initial(g)};{@(g) bb.*g.^(1).*n_initial(g)}];

% SET GROWTH RATES FOR FEM
GR_MC = [{@(v) a*N0*10^Vmin*v.^(-1).*fun(10^Vmin,v)};...
    {@(v) bbb*N0*10^Vmin*v.^(-1).*fun(10^Vmin,v)}];

% tic
% CREATE FEM/PGFEM MATRICES
[A,G1,G2,B,C,~,A_PG,G1_PG,G2_PG,B_PG,C_PG] = ...
    multicomponent_GDE_FE_matrix_creator(g,GR_MC,[],fun,eps);
% toc

% BOUNDARY CONDITIONS FOR FEM
if create_disc == 'y'
    B1 = 10^-Vmin*PHI(1).phi(1,:)';
    B2 = 10^-Vmin*PHI(2).phi(1,:)';
    B1 = [B1;0];
    B2 = [B2;0];
    dB1 = gradient(B1,delta_t);
    dB2 = gradient(B2,delta_t);
    dB1 = [dB1;0];
    dB2 = [dB2;0];
    BC = [{@(ii) B1(ii)},{@(ii) B2(ii)}];
    dBC = [{@(ii) dB1(ii)},{@(ii) dB2(ii)}];
else
    % If discrete is not calculated zero boundary conditions are used.
    % Can be changed to wanted function.
    B1 = [zeros(length(tt),1);0];
    B2 = [zeros(length(tt),1);0];
    dB1 = [zeros(length(tt),1);0];
    dB2 = [zeros(length(tt),1);0];
    BC = [{@(ii) B1(ii)},{@(ii) B2(ii)}];
    dBC = [{@(ii) dB1(ii)},{@(ii) dB2(ii)}];
end

% Calculating the temporal evolutions
[q_FEM_tot,q_FEM,nodes,tt] = ...
    multicomponent_GDE_CrankNicolson(g,initial_dist,delta_t,tmax,...
    A,G1,G2,B,C,[],BC,dBC);


%% MC Sectional method

fun_sec = @(v,w) fun(v,w);
if create_disc == 'y'
[phi_tot_sec,q_tot_sec,PHI_sec,q_sec,v_sec,v_width] = MC_GDE_sectional_method_volume_conc(...
    bins,[log10(min(v_middle)),Vmax],GR_MC,fun,[],[],initial_dist_sec,'volume',tmax,delta_t,[{@(t) 0},{@(t) 0}]);
else
[phi_tot_sec,q_tot_sec,PHI_sec,q_sec,v_sec,v_width] = MC_GDE_sectional_method_volume_conc(...
    bins,[Vmin,Vmax],GR_MC,fun,[],[],initial_dist_sec,'volume',tmax,delta_t,[{@(t) 0},{@(t) 0}]);
end
%% Error computations for the FEM and the sectional method
if create_disc == 'y'
    [error_sec,avg_error_sec] = Error_estimator(v_sec,q_tot_sec,v_middle,10^(-Vmin)*PHI(1).phi+10^(-Vmin)*PHI(2).phi,v_edges,1);
    [error_FEM,avg_error_FEM] = Error_estimator(g,q_FEM_tot,v_middle,10^(-Vmin)*PHI(1).phi+10^(-Vmin)*PHI(2).phi);
end
%% Plotting
if create_disc == 'y'
    % Create figure by using the fig function
    h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13);
    semilogx(v_middle,v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,1)+10^(-Vmin)*PHI(2).phi(:,1)),'k--','LineWidth',1.5)
    hold on
    semilogx(g,g.^(-1).*q_FEM{1}(:,end),'b','LineWidth',1.1)
    semilogx(v_sec,v_sec.^(-1).*q_sec(1).q(:,end),'r','LineWidth',1.1)
    semilogx(v_middle,v_middle.^(-1).*10^(-Vmin).*PHI(1).phi(:,end),'k--','LineWidth',2)
    semilogx(v_middle,v_middle.^(-1).*10^(-Vmin).*PHI(2).phi(:,end),'k--','LineWidth',2)
    semilogx(v_middle,v_middle.^(-1).*10^(-Vmin).*PHI(1).phi(:,1),'k--','LineWidth',1.5)
    semilogx(g,g.^(-1).*q_FEM{2}(:,end),'b','LineWidth',1.1)
    semilogx(v_sec,v_sec.^(-1).*q_sec(2).q(:,end),'r','LineWidth',1.1)
    semilogx(v_middle,v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,end)+10^(-Vmin)*PHI(2).phi(:,end)),'k--','LineWidth',1.5)
    legend('Discrete','FEM','Sectional','location','northwest')
    xlabel('Particle Volume [\mum^3]')
    ylabel('Size distribution [\mum^{-3}cm^{-3}]')
    legend boxoff
    %%
    h2 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13);
    plot(tt,error_FEM,'LineWidth',1.5)
    hold on
    plot(tt,error_sec,'LineWidth',1.5)
    xlabel('Evolution time [s]')
    ylabel('Relative error [%]')
    legend('FEM','Sectional','location','northwest')
    legend boxoff


    % "Banana" plot for the size distribution
    h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13,'border','on');
    [d_grid,t_grid] = meshgrid(v_middle,tt);
    imgrid.d_grid = d_grid;
    imgrid.t_grid = t_grid;
    % labels
    labels.xlab = 'Time [s]';
    labels.ylab = 'Particle Volume [\mum^3]';
    labels.yscale = 'Log';
    labels.ytick = [1e-9 1e-7 1e-6];
    % labels.xtick = [0 10 20 30 40 48];
    labels.xtick = [0 24 48 72 96];
    labels.clab = 'Number distribution [d N/d ln(v_p)]';
    labels.size = 13;
    labels.title = [];
    % Data
    ImagesForScaling = [v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,1)+10^(-Vmin)*PHI(2).phi)];
    fig1.figno = 3;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;

    % Using plotting function, which allows logaritmic axis
    c = PlotParticleDensityEvolution(v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,1)+10^(-Vmin)*PHI(2).phi),imgrid,ImagesForScaling,labels,fig1);
    c.Limits = [0,max(max(v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,1)+10^(-Vmin)*PHI(2).phi)))];
else
    % Create figure by using the fig function
    h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13);
    semilogx(g,g.^(-1).*q_FEM{1}(:,1)+g.^(-1).*q_FEM{2}(:,1),'k--','LineWidth',1.5)
    hold on
    semilogx(g,g.^(-1).*q_FEM{1}(:,end),'b','LineWidth',1.3)
    semilogx(v_sec,v_sec.^(-1).*q_sec(1).q(:,end),'r','LineWidth',1.3)
    semilogx(g,g.^(-1).*q_FEM{2}(:,end),'b','LineWidth',1.3)
    semilogx(v_sec,v_sec.^(-1).*q_sec(2).q(:,end),'r','LineWidth',1.3)
%     semilogx(v_middle,v_middle.^(-1).*(10^(-Vmin)*PHI(1).phi(:,end)+10^(-Vmin)*PHI(2).phi(:,end)),'k--','LineWidth',1.5)
    legend('Initial','FEM','Sectional','location','northwest')
    xlabel('Particle Volume [\mum^3]')
    ylabel('Size distribution [\mum^{-3}cm^{-3}]')
    legend boxoff

    % "Banana" plot for the size distribution
    h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13,'border','on');
    [d_grid,t_grid] = meshgrid(g,tt);
    imgrid.d_grid = d_grid;
    imgrid.t_grid = t_grid;
    % labels
    labels.xlab = 'Time [s]';
    labels.ylab = 'Particle Volume [\mum^3]';
    labels.yscale = 'Log';
    labels.ytick = [1e-9 1e-7 1e-6];
    % labels.xtick = [0 10 20 30 40 48];
    labels.xtick = [0 24 48 72 96];
    labels.clab = 'Number distribution [d N/d ln(v_p)]';
    labels.size = 13;
    labels.title = [];
    % Data
    ImagesForScaling = [g.^(-1).*q_FEM{1}+g.^(-1).*q_FEM{2}];
    fig1.figno = 3;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;

    % Using plotting function, which allows logaritmic axis
    c = PlotParticleDensityEvolution(g.^(-1).*q_FEM{1}+g.^(-1).*q_FEM{2},imgrid,ImagesForScaling,labels,fig1);
    c.Limits = [0,max(max(g.^(-1).*q_FEM{1}+g.^(-1).*q_FEM{2}))];
end


