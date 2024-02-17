%% This file creates different size matrixes for the test problems:

savepath = [pwd,'\GDE_FE_matrices\'];

%% 1 Coagulation matrix creation
% Creates FE matrices for the Coagulation equation

if create(1) == 1
    
    % Initializing cell array
    Coagulation_cell_const = cell(length(sizes),5);
    
    % Initial values for the time evolution 
    beta_const = 8e-6; % Coagulation kernel function value
    
    Vmin = -8; % Minimum size of particle \mu m^3
    Vmax = -2; % Maximum size of particle \mu m^3
    
    % Coagulation kernel function
    fun_const = @(v,w) beta_const;
    
    GR = @(v) 1;
    dGR = @(v) 0;
    
    % Creates coagulation FE matrices corresponding the chosen
    % discretizations
    for kk = 1:length(sizes)
        B = [];
        C = [];
        g = [];
        time = [];
        
        % Node points are spaced logaritmically
        g = logspace(Vmin,Vmax,sizes(kk))';
        g = g(:);
        tic
        [B_const,C_const] = Coagulation_quadrature_matrix_creator(g,fun_const,10^Vmin);
        [~,~,~,~,M] = Condensation_FEM_matrix_creator(...
            sizes(kk),[Vmin Vmax],GR,dGR,0);
        time = toc;
        
        M = cell2mat(M);
        M(1,:) = 0;
        M(:,1) = 0;
        M(1,1) = 1;
        Coagulation_cell_const(kk,:) = [{M},{B_const},{C_const},{g},{time}];
        
        
    end
    clearvars M B_const C_const g
    save([savepath,'Coagulation_FEM_matrices.mat'],'-v7')
end
%% 2 GDE matrix creation
% Creates FE matrices for the test case where both condensation and
% coagulation are affecting the aerosol size distribution (GDE test case 
% with the analytical solution)
if create(2) == 1
    clearvars -except sizes create savepath sizes_cond sizes_discrete
    
    
    % Growth rate and coagulation coefficient
    sigma = 5e-2;
    beta = 2.5e-6; % Coagulation constant
    
    % Size interval
    Vmin = -8;
    Vmax = -2; % Maximum size of particle \mu m^3
    eps = 0.5; % Upwinding factor for Petrov-Galerkin;
    
    % Coagulation kernel function and the growth rate
    fun = @(v,w) beta;
    GR = @(v) sigma*v;
    dGR = @(v) sigma;
    
    % FE matrix creation for the condensation term
    [K_petrov,M_petrov,G,K,M] = Condensation_FEM_matrix_creator(...
        sizes,[Vmin Vmax],GR,dGR,eps);
    
    % FE matrix creation for the coagulation term
    for ii = 1:length(sizes)
        B = [];
        C = [];
        g = [];
        B_petrov = [];
        C_petrov = [];
        time = [];
        g = logspace(Vmin,Vmax,sizes(ii))';
        g = g(:);
        tic
        [B,C,B_petrov,C_petrov]=Coagulation_quadrature_matrix_creator(g,fun,[],eps);
        time = toc;
        Coagulation_cell_GDE(ii,:) = [{B},{C},{B_petrov},{C_petrov},{time}];
        
        K1 = K{ii,1};
        M1 = M{ii,1};
        K1(1,:) = 0;
        K1(:,1) = 0;
        K1(1,1) = 1;
        
        
        M1(1,:) = 0;
        M1(:,1) = 0;
        M1(1,1) = 1;
        
        K{ii,1} = K1;
        M{ii,1} = M1;
        
        K_petrov1 = K_petrov{ii,1};
        M_petrov1 = M_petrov{ii,1};
        K_petrov1(1,:) = 0;
        K_petrov1(:,1) = 0;
        K_petrov1(1,1) = 1;
        
        
        M_petrov1(1,:) = 0;
        M_petrov1(:,1) = 0;
        M_petrov1(1,1) = 1;
        
        K_petrov{ii,1} = K_petrov1;
        M_petrov{ii,1} = M_petrov1;
    end
    save([savepath,'GDE_FEM_matrices.mat'],'Coagulation_cell_GDE','K_petrov','M_petrov','K','M','G','eps','sigma','beta','-v7.3')
end

%% discrete GDE
% FE matrix creation for the test case where the FEM is compared to the
% solution of the discrete GDE
if create(3) == 1
    clearvars -except sizes create savepath sizes_cond sizes_discrete
    
    % Size interval
    Vmin = -9;
    Vmax = -5;
    
    % size of the smallest particle
    vN0 = 1*10^Vmin;
    
    % Value for the size distribution in the smallest bin
    n01 = 1e12;
    % Number of monomers
    N0 = vN0*n01;
    
    % Petrov-Galerkin coefficient
    eps = 0.3;
    
    % Initial conditions
    T = 240; % Temperature K
    rho = 1000; % kg/m^3 Density of particles
    rho = rho*10^-(6*3); % kg/mum^3
    k_b = 1.3806e-23; % m^2*kg*s^-2 K ^-1 Boltzmans constant
    k_b = k_b*10^(2*6); % mum^2...
    b = (3/(4*pi))^(1/6)*sqrt(6*k_b*T/rho)^?;
    
    b = b*10^(-4*1); % unit convertion to cm/s
    
    % Coagulation kernel function
    fun = @(v,w) b*(1./w+1./v).^(1/2).*(v.^(1/3)+w.^(1/3)).^2;
%     sizes_discrete = [20,30,40,50,75,100,125,150,200,250,300,350,400,450,500,600,700,800,900,1000]';
    
    for ii = 1:length(sizes_discrete)
        
        % Smallest particle consist 2 monomers
        g = logspace(log10(2*10^(-9)),Vmax,sizes_discrete(ii))';
        
        % Coagulation matrix creation
        tic
        [B,C,B_petrov,C_petrov] = Coagulation_quadrature_matrix_creator(g,fun,3,eps);
        
        % Condensation matrix creation
        GR = @(v) N0.*10^Vmin*b*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2;
        dGR = @(v) N0.*10^Vmin*b*(-(1/2)*v.^(-2).*(1./(10^Vmin)+1./v).^(-1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2 ...
            +2*(1/3)*v.^(-2/3).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)));
        [K_petrov,M_petrov,~,K,M] = Condensation_FEM_matrix_creator(...
            length(g),g,GR,dGR,eps);
        time = toc;
        
        M = cell2mat(M);
        M_petrov = cell2mat(M_petrov);
        K = cell2mat(K);
        K_petrov = cell2mat(K_petrov);
        
        discrete_GDE_cell(ii,:)=[{g},{B},{C},{B_petrov},{C_petrov},{M},{K},{M_petrov},{K_petrov},{time}];
    end
    % Clearing "unuseful" variables
    clearvars B C B_petrov C_petrov g M M_petrov K_petrov K time

    save([savepath,'discrete_GDE_matrices.mat'],'-v7.3')
    
end
%% 4 Condensation matrix creation
% Creates FE matrices for the condensation equation test case
if create(4) == 1
    clearvars -except sizes create savepath sizes_cond sizes_discrete
%     sizes = [sizes; 1500; 2000; 2500; 3000; 3500; 4000; 4500; 5000];

    % Petrov-Galerkin upwinding factor
    eps = 0.5;
    A = 5e-19*3600; % m^2/h
    tic
    [K_petrov,M_petrov,G,K,M] = Condensation_FEM_matrix_creator(sizes_cond,[-7,-5],@(v) A.*v.^(-1),@(v) -A.*v.^(-2),eps);
    time = toc;
    for ii = 1:length(K_petrov)
        K_petrov1 = K_petrov{ii,1};
        M_petrov1 = M_petrov{ii,1};
        K_petrov1(1,:) = 0;
        K_petrov1(:,1) = 0;
        K_petrov1(1,1) = 1;
        
        
        M_petrov1(1,:) = 0;
        M_petrov1(:,1) = 0;
        M_petrov1(1,1) = 1;
        
        K_petrov{ii,1} = K_petrov1;
        M_petrov{ii,1} = M_petrov1;
        
        K1 = K{ii,1};
        M1 = M{ii,1};
        K1(1,:) = 0;
        K1(:,1) = 0;
        K1(1,1) = 1;
        
        
        M1(1,:) = 0;
        M1(:,1) = 0;
        M1(1,1) = 1;
        
        K{ii,1} = K1;
        M{ii,1} = M1;
    end
    save([savepath,'Condensation_FEM_matrices.mat'],'K','M','K_petrov','M_petrov','G','eps','A','time','-v7.3')
end
