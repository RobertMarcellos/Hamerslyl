function [K_petrov,M_petrov,G,K,M] ...
    = Condensation_FEM_matrix_creator(matrix_sizes,interval,fun,dfun,epsilon)
% [K_petrov,M_petrov,G,time_step_max,K,M,K_tilde,input_time_step_cell] = Condensation_FEM_matrix_creator2(matrix_sizes,interval,fun,dfun,type,epsilon,input_time_step)
%
% Condensation_FEM_matrix_creator calculates Petrov-Galerkin FEM-matrixes
% for condensation equation. Code can be used to calculate multiple
% different size FEM matrixes as once and also logaritmic spacing is
% possible. This version uses different weight matrix for upwinding.
%
% Inputs:
% matrix_sizes  :   vector of wanted matrix sizes (For example [100, 200]
%                   creates FEM matrixes which dimensions are 100x100
%                   and 200x200
% fun           :   Growth rate for aerosol population (I_0)
% dfun          :   Derivative of the growth rate function
% epsilon       :   Upwind factor for Petrov-Galerkin test function
% input_time_step:  Calculates Petrov-Galerkin matrixes otimized for
%                   certain time step which doesn't have to be optimal
%
% Outputs:
% K_petrov      :   Petrov-Galerkin FEM matrix for stiffness matrix
% M_petrov      :   Petrov-Galerkin FEM matrix for mass matrix
% G             :   Cell array for node points
% K             :   Stiffness matrix without Petrov-Galerkin
% M             :   Mass matrix


% Initializing matrixes for outputs
K_petrov = cell(length(matrix_sizes),1);
M_petrov = cell(length(matrix_sizes),1);
K = cell(length(matrix_sizes),1);
M = cell(length(matrix_sizes),1);
G = cell(length(matrix_sizes),1);

for kk = 1:length(matrix_sizes)
    %     g = [];
    %     g = g(:);
    M2 = [];
    K2 = [];
    M1 = [];
    K1 = [];
    
    % Spacing for the node points
    if length(interval) > 2
        g = interval(:);
    else
        g = logspace(interval(1),interval(2),matrix_sizes(kk));
        g = g(:);  
    end
%     end
    

    % Points for the Gaussian quadrature
    ksi = [1/2-sqrt(3)/6 , 1/2+sqrt(3)/6]';
    w = [1/2 1/2]';
    
    
    % Basis functions
    phi1 = 1-ksi;
    phi2 = ksi;
    dphi1 = [-1 -1]';
    dphi2 = [1 1]';
    
    % Petrov-Galerkin test functions
    PGphi1 = epsilon*(6*ksi.^2-6*ksi);
    PGphi2 = epsilon*(-6*ksi.^2+6*ksi);
    
    d_PGphi1 = epsilon*(12*ksi-6);
    d_PGphi2 = epsilon*(-12*ksi+6);
    
    % Matrix from the node points
    H = [1:length(g)-1;2:length(g)]';
    
    [rg,sg]=size(g);
    [rH,sH]=size(H);
    
    % Initialisin the FE matrices
    K2 = zeros(rg,rg); %nollataan harva matriisi K
    M2 = zeros(rg,rg); % Alustetaan massamatriisi
    K1 = zeros(rg,rg); %nollataan harva matriisi K
    M1 = zeros(rg,rg); % Alustetaan massamatriisi
    
    for ii=1:rH
        ind=H(ii,:); % Nodes for the element
        gind=g(ind,:); % Size values for the borders of the node
        
        % Length of the node
        h = gind(2)-gind(1);
        x = h*ksi;
        
        
        % Petrov-Galerkin and Galerkin FEM matrix calculations by using the
        % Gaussian quadrature.
        
        % Mass matrix
        m = [h*sum(w.*phi1.*phi1) h*sum(w.*phi1.*phi2);
            h*sum(w.*phi2.*phi1) h*sum(w.*phi2.*phi2)];
        
        % PG mass matrix
        m2 = [h*sum(w.*phi1.*PGphi1) h*sum(w.*phi1.*PGphi2);
            h*sum(w.*phi2.*PGphi1) h*sum(w.*phi2.*PGphi2)];
        
        % Stiffness matrices
%         k1 = -[h*sum(w.*phi1.*phi1.*dfun(gind(1)+x)) h*sum(w.*phi1.*phi2.*dfun(gind(1)+x));
%             h*sum(w.*phi2.*phi1.*dfun(gind(1)+x)) h*sum(w.*phi2.*phi2.*dfun(gind(1)+x))];
%         
%         k2 = -[sum(w.*dphi1.*phi1.*fun(gind(1)+x)) sum(w.*dphi1.*phi2.*fun(gind(1)+x));
%             sum(w.*dphi2.*phi1.*fun(gind(1)+x)) sum(w.*dphi2.*phi2.*fun(gind(1)+x))];

        k1 = [sum(w.*dphi1.*phi1.*fun(gind(1)+x)) sum(w.*dphi2.*phi1.*fun(gind(1)+x));
            sum(w.*dphi1.*phi2.*fun(gind(1)+x)) sum(w.*dphi2.*phi2.*fun(gind(1)+x))];
%         
        % Petrov-Galerkin stiffness matrices
        
%         k2_petrov = -[sum(w.*dphi1.*fun(gind(1)+x).*PGphi1) sum(w.*dphi1.*fun(gind(1)+x).*PGphi2);
%             sum(w.*dphi2.*fun(gind(1)+x).*PGphi1) sum(w.*dphi2.*fun(gind(1)+x).*PGphi2)];
%         
%         k1_petrov = -[h*sum(w.*phi1.*PGphi1.*dfun(gind(1)+x)) h*sum(w.*phi1.*PGphi2.*dfun(gind(1)+x));
%             h*sum(w.*phi2.*PGphi1.*dfun(gind(1)+x)) h*sum(w.*phi2.*PGphi2.*dfun(gind(1)+x))];
%         
        k1_petrov = [sum(w.*d_PGphi1.*phi1.*fun(gind(1)+x)) sum(w.*d_PGphi2.*phi1.*fun(gind(1)+x));
            sum(w.*d_PGphi1.*phi2.*fun(gind(1)+x)) sum(w.*d_PGphi2.*phi2.*fun(gind(1)+x))];
        
        
        % Local FE matrix
%         k = k1 + k2;
        k = k1;
        
        % Local PGFE matrix
%         k_petrov = k1 + k2 + k1_petrov + k2_petrov;
        k_petrov = k1 + k1_petrov;
        
        % Placing the local matrices into the global matrices
        K1(ind,ind)=K1(ind,ind)+k';
        M1(ind,ind)=M1(ind,ind)+m';
        
        K2(ind,ind)=K2(ind,ind)+k_petrov';
        M2(ind,ind)=M2(ind,ind)+(m+m2)';
        
                  
    end
    
    % Collecting Galerkin and Petrov-Galerkin FEM matrixes into the cell
    % arrays.
    K_petrov{kk,1} = K2;
    M_petrov{kk,1} = M2;
    K{kk,1} = K1;
    M{kk,1} = M1;
    G{kk,1} = g;
end






end



