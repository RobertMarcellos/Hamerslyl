% Calculates the time evolutions and estimation errors for the case of the
% coagulation equation with the constant coagulation kernel.
clearvars -except time_evol_calc
close all
clc

load('Coagulation_FEM_matrices.mat')
saveloc = [pwd,'\Time_evolutions\'];
% Initializing storage variables for the problem
coagulation_time_evolution_const = cell(size(Coagulation_cell_const,1)+1,24);
coagulation_time_evolution_const(1,:) = [{'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'computational_time_FEM'},....
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computational_time_diff'},{'tt'},{'g'},{'d'},{'n_FEM'},{'error_FEM'},{'average_error_FEM'},{'computational_time_diff'},....
    {'n_diff'},{'error_diff'},{'average_error_diff'},{'computational_time_diff'},{'tt'},{'n_teor1'},{'n_teor2'},{'dd'},{'pre_comp_time_sec'}];

%% For accurate analytical solutions

% Initial parameters for the time evolution
V0 = 5e-5;
N0 = 1e4;
V0_const = 2e-5;
N0_const = 1e4;
max_t = 96;
dd = logspace(Vmin,Vmax,3000)';

% Initial distribution
n0 = (N0/V0).*exp(-dd./V0);
n0_const = (N0_const.*dd/V0_const^2).*exp(-dd./V0_const);

% Two choises of time steps
delta_t = [1,0.1]';

% Calculting analytical solutions for the constant coagulation equation
for kk = 1:length(delta_t)
    
    n_teor_const = [];
    tt = [];
    tt = [0:delta_t(kk):max_t]';

for jj = 1:length(tt)
    
    % Analytical solution for the coagulation equation
    if jj == 1
        
        n_teor_const = n0_const(:);
        
    else
        % Constant kernel
        M0 = 2*N0_const/(2+beta_const*N0_const*tt(jj));
        T_const = 1-M0/N0_const;
        n_teor_const(:,jj) = ((N0_const*(1-T_const)^2)/(sqrt(T_const)*V0_const)).*exp(-dd./V0_const).*sinh(dd*sqrt(T_const)./V0_const);
        
    end
    
    % Storing variables into the cell array
    if kk == 1
        coagulation_time_evolution_const{2,21} = n_teor_const;
        coagulation_time_evolution_const{2,9} = tt;
    else
        coagulation_time_evolution_const{2,22} = n_teor_const;
        coagulation_time_evolution_const{2,20} = tt;
    end

    
end
end
coagulation_time_evolution_const{2,23} = dd;



% Looping different discretizations
for kk = 1:size(Coagulation_cell_const,1)
    clearvars -except kk  Coagulation_cell_const coagulation_time_evolution_const beta_const V0_const N0_const Vmin Vmax max_t delta_t saveloc time_evol_calc
    %% Initial values
    % Fetching inialization from pre-made cell array corresponding current
    % test case
    
    % Initial values for coagulation with constant kernel
    M_const = Coagulation_cell_const{kk,1};
    B_const = Coagulation_cell_const{kk,2};
    C_const = Coagulation_cell_const{kk,3};
    g = Coagulation_cell_const{kk,4};
    
    % Initial distribution same as in the previous codes
    num_of_bin = length(g);
    %% Initialization for difference method (linear kernel function)
    d_edges = logspace(Vmin,Vmax,num_of_bin)';
    d_widths = diff(d_edges);
    d = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
    
    % Applying the size splitting operator
    tic
    X = Size_splitting_operator(d);
    coagulation_time_evolution_const{kk+1,24} = toc;
    
    % Initialization for difference method with constant kernel
    n0_diff_const = (N0_const*d./V0_const^2).*exp(-d/V0_const);
    N0_diff_const = n0_diff_const.*d_widths;
    N_diff_const = N0_diff_const(:);
    n_diff_const = n0_diff_const(:);
    
    coagulation_time_evolution_const{kk+1,10} = g;
    coagulation_time_evolution_const{kk+1,11} = d;
    
%% Time evolution calculations here
    for ii = 1:length(delta_t)
        tt = [0:delta_t(ii):max_t]';
        
        %% Difference method time evolution and error
        N_diff_const = [];
        n_diff_const = [];
        N_diff_const = N0_diff_const(:);
        n_diff_const = n0_diff_const(:);
        
        
        tic
        for jj = 1:length(tt)-1
            F = [];
            R = [];
            % Formation matrix F creation
            F = zeros(length(N0_diff_const));
            for h = 1:length(N0_diff_const)
                
                F(h,:) = N_diff_const(:,jj)'*X{h};
                
                
            end
            F = (beta_const/2)*F;
            R = beta_const*ones(length(n0_diff_const)).*N_diff_const(:,jj);
            
            % Crank-Nicolson for difference method
            N_diff_const(:,jj+1) = (eye(size(F))-delta_t(ii)/2*(F-R))\(N_diff_const(:,jj)+delta_t(ii)/2*(F-R)*N_diff_const(:,jj));
            n_diff_const(:,jj+1) = N_diff_const(:,jj+1)./d_widths;
            
            
        end
        time_diff_const = toc;
        
        % Error calculations
        dd = coagulation_time_evolution_const{2,23};
        if ii == 1
            n_teor = coagulation_time_evolution_const{2,21};
        elseif ii == 2
            n_teor = coagulation_time_evolution_const{2,22};
        end
        
        % Error calculation with Error_esimator function
        [coagulation_time_evolution_const{kk+1,(ii-1)*11+6},coagulation_time_evolution_const{kk+1,(ii-1)*11+7}] = Error_estimator(d,n_diff_const,dd,n_teor,d_edges,1);
        
        % Assigning obtained values into the matrix
        coagulation_time_evolution_const{kk+1,(ii-1)*11+5} = n_diff_const;
        coagulation_time_evolution_const{kk+1,(ii-1)*11+8} = time_diff_const;
        
        
        %% Finite element calculations and errors
         %% FEM initializations and time evolutions
        n0_FEM_const = (N0_const*g./(V0_const^2)).*exp(-g/V0_const);
        n_FEM_const = [];
        n_FEM_const = n0_FEM_const(:);
        [ n_FEM_const,~,time_FEM_const ] = CrankNicolsonGDE(M_const,zeros(size(M_const)),B_const,C_const,n_FEM_const,[0,max_t],delta_t(ii),[],[]);
        
        % Error estimation with Error_estimator
        [coagulation_time_evolution_const{kk+1,(ii-1)*11+2},coagulation_time_evolution_const{kk+1,(ii-1)*11+3}] = Error_estimator(g,n_FEM_const,dd,n_teor);
  
        % Assigning values into the cell array
        coagulation_time_evolution_const{kk+1,(ii-1)*11+1} = n_FEM_const;
        coagulation_time_evolution_const{kk+1,(ii-1)*11+4} = time_FEM_const;     
        
    end
    
    
end


save([saveloc,'\coagulation_time_evolutions.mat'],'coagulation_time_evolution_const','-v7.3')




