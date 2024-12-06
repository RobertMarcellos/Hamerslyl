% This code calculates the discrete GDE and FE and sectional approximation
% with Brownian coagulation kernel
clearvars -except time_evol_calc create_discrete
close all
clc
load('discrete_GDE_matrices.mat')
saveloc = [pwd,'\Time_evolutions\'];

% Initialization for the cell array
discrete_GDE_evolutions_numerical = cell(size(discrete_GDE_cell,1)+1,16);
discrete_GDE_evolutions_numerical(1,:) = [{'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'FEM_comp_time'},...
    {'n_PGFEM'},{'error_PGFEM'},{'average_error_PGFEM'},{'PGFEM_comp_time'},...
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'diff_comp_time'},...
    {'tt'},{'nodes'},{'bins'},{'pre_comp_time_sec'}];

% Initial size distribution
% Creating initial distribution
vv = 10^Vmin:10^Vmin:10^Vmax;
vv = vv(:);
sigma1 = 3e-8;
mu1 = 40*10^Vmin;
N1 = 9e6;
f1 = N1*exp(-0.5*((vv-mu1)./sigma1).^2);

sigma2 = 2e-7;
mu2 = 3e-7;
N2 = 9e6;

f2 = N2*exp(-0.5*((vv-mu2)./sigma2).^2);
n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2)+N2*exp(-0.5*((v-mu2)./sigma2).^2);

% Time evolution parameters
tmax = 0.5;
delta_t = 0.01;
tt = [0:delta_t:tmax]';
discrete_GDE_evolutions_numerical{2,13} = tt;
%% Discrete GDE
% Creation for the discrete GDE
if create_discrete == 'y'
    tic
    [N_disc,n_disc,d_disc,dd_disc,tt_disc] = discrete_GDE_solver_ver2(vv(2),10^Vmax,fun, n_initial,@(t) 0, delta_t, tmax, 'crank',vN0,n01);
    time = toc;
    save([saveloc,'Discrete_evolution_GDE_test_case.mat'], 'N_disc','n_disc','d_disc','dd_disc','tt_disc','vv','time')
else
    load('Discrete_evolution_GDE_test_case.mat')
end

for kk = 1:size(discrete_GDE_cell,1)
    
    % Loading variable for time evolution calculations
    g = discrete_GDE_cell{kk,1};
    B = discrete_GDE_cell{kk,2};
    C = discrete_GDE_cell{kk,3};
    M = discrete_GDE_cell{kk,6};
    K = discrete_GDE_cell{kk,7};
    
    B_petrov = discrete_GDE_cell{kk,4};
    C_petrov = discrete_GDE_cell{kk,5};
    M_petrov = discrete_GDE_cell{kk,8};
    K_petrov = discrete_GDE_cell{kk,9};
    
    
    % FE/PGFE approximations for the time evolutions
    [n_FEM,~,FEM_comp_time] = CrankNicolsonGDE( M,K,B,C,n_initial(g),[0,tmax],delta_t,[],[]);
    % Error of the FE approximation
    [discrete_GDE_evolutions_numerical{kk+1,2},discrete_GDE_evolutions_numerical{kk+1,3}] = Error_estimator(g,n_FEM,d_disc,n_disc);
    discrete_GDE_evolutions_numerical{kk+1,1} = n_FEM;
    discrete_GDE_evolutions_numerical{kk+1,4} = FEM_comp_time;
    
    [n_PGFEM,~,PGFEM_comp_time] = CrankNicolsonGDE( M_petrov,K_petrov,B_petrov,C_petrov,n_initial(g),[0,tmax],delta_t,[],[]);
    % Error of the PGFE approximation
    [discrete_GDE_evolutions_numerical{kk+1,6},discrete_GDE_evolutions_numerical{kk+1,7}] = Error_estimator(g,n_PGFEM,d_disc,n_disc);
    discrete_GDE_evolutions_numerical{kk+1,5} = n_PGFEM;
    discrete_GDE_evolutions_numerical{kk+1,8} = PGFEM_comp_time;
    discrete_GDE_evolutions_numerical{kk+1,14} = g;
    
    % Sectional method approximations
    clearvars d_edges d_widths d R X nd bb a I it g_true Beta Beta_form AA gg AAtmp n_diff N_diff F
    d_edges = logspace(log10(min(d_disc)),Vmax,length(g))';
    d_widths = diff(d_edges);
    d = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
    R = zeros(length(d));
    
    tic
    % Size splitting operator for the coagulation
    X = Size_splitting_operator(d);
    
    % Initializations for the sectional method matrices
    nd = length(d);
    bb = -1./diff(d);
    a = -(d_widths(2:end)./d_widths(1:end-1)).*bb;
    bb = [bb(1); bb];
    GR = @(v) N0.*10^Vmin*b*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2;
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
    discrete_GDE_evolutions_numerical{kk+1,16} = toc;
    
    
    % Initial distribution
    n_diff(:,1) = n_initial(d);
    N_diff(:,1) = n_diff.*d_widths;
    
    
    tic
    % Time evolution for the sectional method
    for ii = 1:length(tt)-1
        
        % Formation matrix F creation
        F = zeros(length(d));
        for j = 1:length(d)
            
            F(j,:) = N_diff(:,ii)'*(0.5*Beta_form.*X{j});
            
            
        end
        R = N_diff(:,ii).*Beta;
        % Crank-Nicolson for difference method
        N_diff(:,ii+1) = (eye(size(F))-delta_t/2*(F-R+AA))\(N_diff(:,ii)+delta_t/2*(F-R+AA)*N_diff(:,ii));
        n_diff(:,ii+1) = N_diff(:,ii+1)./d_widths;
    end
    diff_comp_time = toc;
    discrete_GDE_evolutions_numerical{kk+1,9} = n_diff;
    discrete_GDE_evolutions_numerical{kk+1,12} = diff_comp_time;
    [discrete_GDE_evolutions_numerical{kk+1,10},discrete_GDE_evolutions_numerical{kk+1,11}] = Error_estimator(d,n_diff,d_disc,n_disc,d_edges,1);
    discrete_GDE_evolutions_numerical{kk+1,15} = d;
    
    
    
end
save([saveloc,'discrete_GDE_test_case.mat'],'discrete_GDE_evolutions_numerical')