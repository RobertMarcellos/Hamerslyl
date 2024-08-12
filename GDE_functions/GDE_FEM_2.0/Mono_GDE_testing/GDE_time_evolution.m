% This code calculates GDE time evolutions with constant coagulation kernel
clearvars -except time_evol_calc 
close all
clc
load('GDE_FEM_matrices.mat')
saveloc = [pwd,'\Time_evolutions\'];

% Initialization for the cell array
GDE_time_evolutions = cell(size(G,1)+1,31);
GDE_time_evolutions(1,:) = [{'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'computation_time_FEM'},...
    {'n_PGFEM'},{'error_PGFEM'},{'average_error_PGFEM'},{'computation_time_PGFEM'},...
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computation_time_diff'},...
    {'tt1'},{'nodes'},{'bins'},...
    {'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'computation_time_FEM'},...
    {'n_PGFEM'},{'error_PGFEM'},{'average_error_PGFEM'},{'computation_time_PGFEM'},...
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computation_time_diff'},...
    {'tt2'},{'n_teor1'},{'n_teor2'},{'dd'}];

% Initialization for the analytical solution
N0 = 1e4; % Number of particles
V0 = 2e-7; % Mean particle volume
Vmin = log10(min(G{1})); % Minimum size of particle \mu m^3
Vmax = log10(max(G{1})); % Maximum size of particle \mu m^3
max_t = 96;

% Decisions for the time steps
delta_t = [1 0.1]';

% Loop for the time steps
for ii = 1:length(delta_t)
    
    % Analytical time evolutions
    dd = logspace(Vmin,Vmax,3000)';
    n0 = (N0/V0)*exp(-dd./V0);
    tt = [0:delta_t(ii):max_t]';
    
    for kk = 1:length(tt)
        if kk == 1
            n_teor = (N0*dd/V0^2).*exp(-dd/V0);
            
        else
            M0 = 2*N0/(2+beta*N0*tt(kk));
            M1 = N0*V0*exp(sigma*tt(kk));
            n_teor(:,kk) = M0^2/(M1*sqrt(1-M0/N0))*exp(-N0*dd/M1).*sinh(sqrt(1-M0/N0)*N0*dd/M1);
            
        end
        
        
        
        
    end
    
    if ii == 1
        GDE_time_evolutions{2,29} = n_teor;
        GDE_time_evolutions{2,13} = tt;
        GDE_time_evolutions{2,31} = dd;
    elseif ii == 2
        GDE_time_evolutions{2,30} = n_teor;
        GDE_time_evolutions{2,28} = tt;
    end
    
    for kk = 1:length(G)
        
        %% Initializations for the FE approximations
        g = G{kk};
        n0_FEM = (N0*g/V0^2).*exp(-g/V0);
        K_loc = K{kk};
        M_loc = M{kk};
        B = Coagulation_cell_GDE{kk,1};
        C = Coagulation_cell_GDE{kk,2};
        
        K_petrov_loc = K_petrov{kk};
        M_petrov_loc = M_petrov{kk};
        B_petrov = Coagulation_cell_GDE{kk,3};
        C_petrov = Coagulation_cell_GDE{kk,4};
        
        % FE approximation
        [n_FEM,~,FEM_comp_time] = CrankNicolsonGDE(M_loc,K_loc,B,C,n0_FEM,[0 max_t]',delta_t(ii),@(t) 0,@(t) 0);
        [GDE_time_evolutions{kk+1,2+(ii-1)*15},GDE_time_evolutions{kk+1,3+(ii-1)*15}] = Error_estimator(g,n_FEM,dd,n_teor);
        GDE_time_evolutions{kk+1,1+(ii-1)*15} = n_FEM;
        GDE_time_evolutions{kk+1,4+(ii-1)*15} = FEM_comp_time;
        
        
        % PGFE approximation
        [n_PGFEM,~,PGFEM_comp_time] = CrankNicolsonGDE(M_petrov_loc,K_petrov_loc,B_petrov,C_petrov,n0_FEM,[0 max_t]',delta_t(ii),@(t) 0,@(t) 0);
        [GDE_time_evolutions{kk+1,6+(ii-1)*15},GDE_time_evolutions{kk+1,7+(ii-1)*15}] = Error_estimator(g,n_PGFEM,dd,n_teor);
        GDE_time_evolutions{kk+1,5+(ii-1)*15} = n_PGFEM;
        GDE_time_evolutions{kk+1,8+(ii-1)*15} = PGFEM_comp_time;
        
        GDE_time_evolutions{kk+1,14} = g;
        
        %% Difference approximation
        
        % Initializations for the difference method
        num_of_bin = length(g);
        
        d_edges = logspace(Vmin,Vmax,num_of_bin)';
        d_widths = diff(d_edges);
        d = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
        nd = length(d);
        b = -1./diff(d);
        a = -(d_widths(2:end)./d_widths(1:end-1)).*b;
        b = [b(1); b];
        I = zeros(nd,1);
        I = sigma*(d).^(1);
        it = 1;
        g_true = [];
        g_true(:,it) = I;
        AA = diag(b.*g_true(:,it));
        gg = g_true(:,it);
        AAtmp = diag(a.*gg(1:end-1));
        AA(2:end,1:end-1) = AA(2:end,1:end-1) + AAtmp;
        
        X = Size_splitting_operator(d);
        n0_diff = (N0*d/V0^2).*exp(-d/V0);
        N0_diff =n0_diff.*d_widths;
        N_diff = N0_diff(:);
        n_diff = n0_diff(:);
        
        tic
        for jj = 1:length(tt)-1
            
            % Difference method calculations
            
            F = [];
            R = [];
            %     beta_F = [];
            
            % Formation matrix F creation
            F = zeros(length(N0_diff));
            for h = 1:length(N0_diff)
                
                F(h,:) = N_diff(:,jj)'*X{h};
                
                
            end
            F = (beta/2)*F;
            R = beta*ones(length(n0_diff)).*N_diff(:,jj);
            
            % Crank-Nicolson for difference method
            N_diff(:,jj+1) = (eye(size(F))-delta_t(ii)/2*(F-R+AA))\(N_diff(:,jj)+delta_t(ii)/2*(F-R+AA)*N_diff(:,jj));
            n_diff(:,jj+1) = N_diff(:,jj+1)./d_widths;
            
            
        end
        time_diff = toc;
        % Error calculation
        [GDE_time_evolutions{kk+1,10+(ii-1)*15},GDE_time_evolutions{kk+1,11+(ii-1)*15}] = Error_estimator(d,n_diff,dd,n_teor,d_edges,1);
        
        GDE_time_evolutions{kk+1,9+(ii-1)*15} = n_diff;
        GDE_time_evolutions{kk+1,12+(ii-1)*15} = time_diff;
        GDE_time_evolutions{kk+1,15} = d;
        
    end
    
    
    
end
save([saveloc,'\GDE_time_evolutions.mat'],'GDE_time_evolutions','-v7.3')