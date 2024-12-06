% This function calculates the time evolution for the coagulation equation
% This code is used to calculate time evolution of condensation equation
% using FEM and difference method. Solution is compared to analytical
% solution in a special case.
clearvars -except time_evol_calc
close all
clc

% Loading premade matrixes
load('Condensation_FEM_matrices.mat')
saveloc = [pwd,'\Time_evolutions\'];
cond_time_evolution = cell(size(G,1)+1,37);
cond_time_evolution(1,:) = [{'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'computational_time_FEM'},...
    {'n_PGFEM'},{'error_PGFEM'},{'average_error_PGFEM'},{'computational_time_FEM'},...
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computational_time_diff'},...
    {'time_steps1'},{'elements'},{'bins'},...
    {'n_PGFEM'},{'error_PGFEM'},{'average_error_PGFEM'},{'computational_time_FEM'},...
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computational_time_diff'},...
    {'time_steps2'},{'n_teor_t1'},{'n_teor_t2'},{'dd'},...
    {'n_FEM_ie'},{'error_FEM_ie'},{'average_error_FEM_ie'},{'computational_time_FEM_ie'},...
    {'n_PGFEM_ie'},{'error_PGFEM_ie'},{'average_error_PGFEM_ie'},{'computational_time_FEM_ie'},{'error_sec_mod2_t1'},{'error_sec_mod2_t2'}];

% Initializations for problem
mean_d = 0.3e-6;
mu = 1.2;
N0 = 180;
Vmin = log10(min(G{1}));
Vmax = log10(max(G{1}));
eval_time = 96;

% Time step for time evolutions
delta_t = [1,0.1]';

for ii = 1:length(delta_t)
    
    tt = [0:delta_t(ii):eval_time]';
    
    % Analytical solutions for the test case with same time steps
    dd = logspace(Vmin,Vmax,5000)';
    n_teor = [];
    index2 = [];
    for jj = 1:length(tt)
        
        index2 = find(dd <= sqrt(A*2*tt(jj)));
        
        d2 = dd;
        d2(index2) = 0;
        
        n_teor(:,jj) = (d2./sqrt(d2.^2-2*A*tt(jj))).*(N0./(sqrt(2*pi).*sqrt(d2.^2-2*A*tt(jj))*log(mu))).*exp(-log(sqrt(d2.^2-2*A*tt(jj))./mean_d).^2/(2*log(mu)^2));
        
    end
    
    % Storing variables into the cell array
    if ii == 1
        cond_time_evolution{2,25} = n_teor;
        cond_time_evolution{2,13} = tt;
    elseif ii == 2
        cond_time_evolution{2,26} = n_teor;
        cond_time_evolution{2,24} = tt;
    end
    cond_time_evolution{2,27} = dd;
    
    
    
    for kk = 1:length(K_petrov)
        
        % Initializing FE/PGFEM approximation
        M_loc = [];
        K_loc = [];
        
        M_loc = M{kk,1};
        K_loc = K{kk,1};
        
        
        K_petrov_loc = [];
        M_petrov_loc = [];
        g = [];
        

        K_petrov_loc = K_petrov{kk,1};
        M_petrov_loc = M_petrov{kk,1};
        g = G{kk,1};
        
        % Initial size distribution
        n0_FEM = (N0./(sqrt(2*pi).*g*log(mu))).*exp(-log(g./mean_d).^2/(2*log(mu)^2));
        
        % Time evolution calculations for FE/PGFE approximation
        [n_PGFEM,~,PGFEM_comp_time] = TimeEvolutionGDE(M_petrov_loc,K_petrov_loc,[],[],n0_FEM,[0,eval_time]',delta_t(ii),@(t) 0,@(t) 0,[],'crank');
        [cond_time_evolution{kk+1,6+(ii-1)*11},cond_time_evolution{kk+1,7+(ii-1)*11}] = Error_estimator(g,n_PGFEM,dd,n_teor);
        cond_time_evolution{kk+1,5+(ii-1)*11} = n_PGFEM;
        cond_time_evolution{kk+1,8+(ii-1)*11} = PGFEM_comp_time;
        
        % If time step is 1h
        
        % Error computations
        if ii == 1
            
            [n_FEM,~,FEM_comp_time] = TimeEvolutionGDE(M_loc,K_loc,[],[],n0_FEM,[0,eval_time]',delta_t(ii),@(t) 0,@(t) 0,[],'crank');
            [cond_time_evolution{kk+1,2},cond_time_evolution{kk+1,3}] = Error_estimator(g,n_FEM,dd,n_teor);
            cond_time_evolution{kk+1,1} = n_FEM;
            cond_time_evolution{kk+1,4} = FEM_comp_time;
            
            [n_FEM_ie,~,FEM_comp_time_ie] = TimeEvolutionGDE(M_loc,K_loc,[],[],n0_FEM,[0,eval_time]',delta_t(ii),@(t) 0,@(t) 0,[],'impl.euler');
            [cond_time_evolution{kk+1,29},cond_time_evolution{kk+1,30}] = Error_estimator(g,n_FEM_ie,dd,n_teor);
            cond_time_evolution{kk+1,28} = n_FEM_ie;
            cond_time_evolution{kk+1,31} = FEM_comp_time_ie;
            
            [n_PGFEM_ie,~,PGFEM_comp_time_ie] = TimeEvolutionGDE(M_petrov_loc,K_petrov_loc,[],[],n0_FEM,[0,eval_time]',delta_t(ii),@(t) 0,@(t) 0,[],'impl.euler');
            [cond_time_evolution{kk+1,33},cond_time_evolution{kk+1,34}] = Error_estimator(g,n_PGFEM_ie,dd,n_teor);
            cond_time_evolution{kk+1,32} = n_PGFEM_ie;
            cond_time_evolution{kk+1,35} = PGFEM_comp_time_ie;
        end
        cond_time_evolution{kk+1,14} = g;
        
        
        % Initializing sectional method solution
        d_edges = g;
        d_widths = diff(d_edges);
        d = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
        nd = length(d);
        b = -1./diff(d);
        a = -(d_widths(2:end)./d_widths(1:end-1)).*b;
        b = [b(1); b];
        g_true = zeros(nd,1);
        g_true_tmp = A*(d).^(-1);
        %     g_true_tmp(1) = 1e-4;
        it = 1;
        g_true(:,it) = g_true_tmp;
        AA = diag(b.*g_true(:,it));
        gg = g_true(:,it);
        AAtmp = diag(a.*gg(1:end-1));
        AA(2:end,1:end-1) = AA(2:end,1:end-1) + AAtmp;
        
        % Initial distribution
        n0_diff = (N0./(sqrt(2*pi).*d*log(mu))).*exp(-log(d./mean_d).^2/(2*log(mu)^2));
        N0_diff = d_widths.*n0_diff;
        
        % Sectional method time evolution
        [N_diff,~,diff_comp_time] = TimeEvolutionGDE(eye(size(AA)),AA,[],[],N0_diff,[0,eval_time]',delta_t(ii),@(t) 0,@(t) 0,[],'crank');
        n_diff = N_diff./d_widths;
        % Error calculation
        [cond_time_evolution{kk+1,10+(ii-1)*11},cond_time_evolution{kk+1,11+(ii-1)*11}] = Error_estimator(d,n_diff,dd,n_teor,d_edges,1);
        [cond_time_evolution{kk+1,36+(ii-1)}] = Error_estimator(d,n_diff,dd,n_teor,d_edges,2);
        cond_time_evolution{kk+1,9+(ii-1)*11} = n_diff;
        cond_time_evolution{kk+1,12+(ii-1)*11} = diff_comp_time;
        cond_time_evolution{kk+1,15} = d;

        
        
    end
    
    
    
    
end
save([saveloc,'\condensation_time_evolutions.mat'],'cond_time_evolution','-v7.3')



