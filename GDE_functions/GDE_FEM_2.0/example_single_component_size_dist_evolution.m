% This code calculates an example temporal evolution with chosen parameter
% for monodisperse General Dynamic Equation of Aerosols (GDE), which is
% described in detail in Seinfeld, J. H., & Pandis, S. N. (2016).
% Atmospheric chemistry and physics: From air pollution to climate change,
% 3rd edition.
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

% CREATE DISCRETE GDE TEMPORAL EVOLUTION
create_disc = 'y'; % 'n'/'y' (no/yes)

% PLOT INITIAL SIZE DISTRIBUTION
plot_init = 'y'; % 'n'/'y' (no/yes)

% NUMBER OF BINS AND ELEMENTS
bins = 100;
elements = 100;

% PARAMETERS FOR THE SIZE INTERVAL
Vmin = -9;
Vmax = -6; % The size interval is chosen to be narrow because of discrete GDE

% PETROV-GALERKIN PARAMETER
eps = 0.3;

% Discretization for the discrete GDE
vv = 10^Vmin:10^Vmin:10^Vmax;
vv = vv(:);

% Changing the value of the number distribution to value of number
% concentration (for monomers)
vN0 = vv(1);
n01 = 5e10;
N0 = vN0*n01;

% CREATING INITIAL DISTRIBUTION
sigma1 = 2e-8;
mu1 = 50*10^Vmin;
N1 = 5e6;
f1 = N1*exp(-0.5*((vv-mu1)./sigma1).^2);
n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2);
n0 = n_initial(vv);

% TEMPORAL DISCRETIZATION PARAMETER
tmax = 96;
delta_t = 0.5;
tt = [0:delta_t:tmax]';

if plot_init == 'y'
    figure
    semilogx(vv,n0,'k','LineWidth',1.7)
    xlabel('Particle volume [\mum^3]')
    ylabel('Size distribution [[\mum^{-3}cm^{-3}]')
    drawnow
end

% CREATE COAGULATION KERNEL
% This can be changed. Here Brownian coagulation kernel is applied

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
if create_disc == 'y'
    [N_disc,n_disc,d_disc,dd_disc,tt_disc] = discrete_GDE_solver_final(vv(2),10^Vmax,fun, n_initial,@(t) 0,[], delta_t, tmax, 'crank',vN0,n01,'GDE');
end
%% FEM matrix creations for the problem
if create_disc == 'y'
    g = logspace(log10(min(d_disc)),Vmax,elements)';
else
    g = logspace(Vmin,Vmax,elements)';
end
g = g(:);

% CREATE COAGULATION FEM MATRICES
% Inputs: 1. node points, 2. coagulation kernel, 3. number of quadrature
% points (3,5), Petrov-Galerkin parameter
[B,C,B_petrov,C_petrov] = Coagulation_quadrature_matrix_creator2(g,fun,3,eps);

% CREATE CONDENSATION FEM MATRICES
% Function for condensation can be changed!
GR = @(v) 3600*10^(-9)*N0.*10^Vmin*b*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2;
% Derivative of GR
dGR = @(v) 3600*10^(-9)*N0.*10^Vmin*b*(-(1/2)*v.^(-2).*(1./(10^Vmin)+1./v).^(-1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2 ...
    +2*(1/3)*v.^(-2/3).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)));

if exist('d_disc','var')
    [K_petrov,M_petrov,~,K,M] = Condensation_FEM_matrix_creator(...
        elements,[log10(min(d_disc)),Vmax],GR,dGR,eps);
else
    [K_petrov,M_petrov,~,K,M] = Condensation_FEM_matrix_creator(...
        elements,[Vmin,Vmax],GR,dGR,eps);
end
M = cell2mat(M);
M_petrov = cell2mat(M_petrov);
K = cell2mat(K);
K_petrov = cell2mat(K_petrov);

%% FE time evolutions
n_FEM(:,1) = n_initial(g);
n_PGFEM(:,1) = n_initial(g);
if exist('n_disc','var')
    derivative = [0;diff(n_disc(1,:))'];
    [n_FEM,tt1] = CrankNicolsonGDE_new( M,K,B,C,n_FEM,[0,tmax],delta_t,n_disc(1,:)',derivative);
else
    % With zero boundary condition
    [n_FEM,tt1] = CrankNicolsonGDE_new( M,K,B,C,n_FEM,[0,tmax],delta_t,zeros(length(tt),1),zeros(length(tt),1));
end

%% Sectional method for the GDE

if create_disc == 'y'
    d_edges = logspace(log10(min(d_disc)),Vmax,bins)';
else
    d_edges = logspace(Vmin,Vmax,bins)';
end
d_widths = diff(d_edges);
d = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
R = zeros(length(d));

% Creating size splitting operator for the sectional method
X = Size_splitting_operator(d);

nd = length(d);


bb = -1./diff(d);
a = -(d_widths(2:end)./d_widths(1:end-1)).*bb;
bb = [bb(1); bb];

% GR = @(v) 3600*10^(-9)*N0.*10^Vmin*b*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2;
I = zeros(nd,1);
I = GR(d);

it = 1;
g_true = [];
g_true(:,it) = I;

AA = diag(bb.*g_true(:,it));
gg = g_true(:,it);
AAtmp = diag(a.*gg(1:end-1));
AA(2:end,1:end-1) = AA(2:end,1:end-1) + AAtmp;

for j = 1:length(d)

    Beta(:,j) = fun(d,d(j));

end
Beta_form = Beta;

% Initial distribution
n_diff(:,1) = n_initial(d);
N_diff(:,1) = n_diff.*d_widths;
tic
for kk = 1:length(tt)-1

    % Formation matrix F creation
    F = zeros(length(d));
    for j = 1:length(d)

        F(j,:) = N_diff(:,kk)'*(0.5*Beta_form.*X{j});


    end
    R = N_diff(:,kk).*Beta;

    % Crank-Nicolson for difference method
    N_diff(:,kk+1) = (eye(size(F))-delta_t/2*(F-R+AA))\(N_diff(:,kk)+delta_t/2*(F-R+AA)*N_diff(:,kk));
    n_diff(:,kk+1) = N_diff(:,kk+1)./d_widths;
end

%% Error computations for the FEM and the sectional method
if create_disc == 'y'
    [error_sec,avg_error_sec] = Error_estimator(d,n_diff,d_disc,n_disc,d_edges,1);
    [error_FEM,avg_error_FEM] = Error_estimator(g,n_FEM,d_disc,n_disc);
end
%% Plotting

% Create figure by using the fig function
h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13);

if create_disc == 'y'
    semilogx(d_disc,n_disc(:,1),'k--','LineWidth',1.5)
    hold on
    semilogx(g,n_FEM(:,end),'b','LineWidth',1.1)
    semilogx(d,n_diff(:,end),'r','LineWidth',1.1)
    semilogx(d_disc,n_disc(:,end),'k--','LineWidth',1.5)
    legend('Discrete','FEM','Sectional','location','northwest')
    xlabel('Particle Volume [\mum^3]')
    ylabel('Size distribution [\mum^{-3}cm^{-3}]')
    legend boxoff

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
    [d_grid,t_grid] = meshgrid(d_disc,tt_disc);
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
    ImagesForScaling = [n_disc];
    fig1.figno = 5;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;
    c = PlotParticleDensityEvolution(n_disc,imgrid,ImagesForScaling,labels,fig1);
    c.Limits = [0,max(max(n_disc))];
else
    semilogx(g,n_FEM(:,1),'k--','LineWidth',1.5)
    hold on
    semilogx(g,n_FEM(:,end),'b','LineWidth',1.1)
    semilogx(d,n_diff(:,end),'r','LineWidth',1.1)
    legend('Initial','FEM','Sectional','location','northwest')
    xlabel('Particle Volume [\mum^3]')
    ylabel('Size distribution [\mum^{-3}cm^{-3}]')
    legend boxoff

    % "Banana" plot for the size distribution
    h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',13,'border','on');
    [d_grid,t_grid] = meshgrid(g,tt1);
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
    ImagesForScaling = [n_FEM];
    fig1.figno = 5;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;
    c = PlotParticleDensityEvolution(n_FEM,imgrid,ImagesForScaling,labels,fig1);
    c.Limits = [0,max(max(n_FEM))];
end




