% This code calculates the discrete tive evolutions with chosen parameters.
clearvars -except discretizations time_steps create_disc_coag create_disc_GDE time_steps_disc
close all
clc

Vmin = -9;
Vmax = -5;

% size of the smallest particle
vN0 = 1*10^Vmin;

% Value for the size distribution in the smallest bin
n01 = 1e16;
% Number of monomers
N0 = vN0*n01;

% Initial conditions
T = 280; % Temperature K
rho = 2*10^(-15); % kg/mum^3
k_b = 1.3806e-23; % m^2*kg*s^-2 K ^-1 Boltzmans constant
k_b = k_b*10^(2*6); % mum^2...
vakio = (3/(4*pi))^(1/6)*sqrt(6*k_b*T/rho);
muunnin = 12; % unit convertion to cm/s

% Brownian coagulation kernel function
fun = @(v,w) 3600*10^(-muunnin)*(vakio*(1./w+1./v).^(1/2).*(v.^(1/3)+w.^(1/3)).^2);

% Initial size distribution
% Creating initial distribution
vv = 10^Vmin:10^Vmin:10^Vmax;
vv = vv(:);
sigma1 = 2e-8;
mu1 = 60*10^Vmin;
N1 = 3e12;
f1 = N1*exp(-0.5*((vv-mu1)./sigma1).^2);

sigma2 = 0.6e-7;
mu2 = 1.5e-7;
N2 = 0.8e12;

f2 = N2*exp(-0.5*((vv-mu2)./sigma2).^2);
n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2)+N2*exp(-0.5*((v-mu2)./sigma2).^2);

% Time evolution parameters
tmax = 2;
delta_t = time_step_disc(1);
tt = [0:delta_t:tmax]';


%% SINGLE COMPONENT
% Single component evolution for comparison and to confirm that 
% multivolume GDE works
tic 
[N_disc,n_disc,d_disc,dd_disc,tt_disc] = discrete_GDE_solver_final(vv(2),10^Vmax,fun, n_initial,[],[], delta_t, tmax+delta_t, 'crank',vN0,n01,'coag');
discrete_GDE_evolution_time = toc;

%% MULTICOMPONENT
% Fractions for the monomers
a = 0.85;
bbb = 0.15;

% Fractions for the initial size distribution

% Dividing initial size distribution to two components
x1 = 1e-9;
x2 = 2e-7;
y1 = 0.8;
y2 = 0.2;
k = (y2-y1)/(x2-x1);
b = y1-k*x1;
aa = k*vv(2:end)+b;
aa(aa < y2) = y2;
bb = 1-aa;

% Multivolume discrete GDE
tic
[N,PHI,v_middle,v_edges,tt] ...
    = discrete_multivolume_GDE(10^Vmin,[10^Vmin,10^Vmax],@(vv) n_initial(vv),'distribution',[{aa},{bb}],...
    [a*vN0*n01,bbb*vN0*n01],fun,delta_t,tmax+delta_t,'coag');
MC_discrete_GDE_evolution_time = toc;

% Save discrete time evolutions
save([pwd,'/Time_evolutions/MC_discrete_GDE_coagulation_only.mat'])



