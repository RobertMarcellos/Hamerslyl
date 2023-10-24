function [ n_FEM,tt,evolution_time ] = CrankNicolsonGDE( M,K,B,C,n_init,t_interval,h,C_in,dC_in,source )
%CRANKNICOLSONGDE Solves FEM matrix equation for whole GDE by using the
%Crank-Nicolson method.

%   INPUTS:
%   M : mass matrix (matrix in the left hand side of FEM equation)
%   K : Stiffness matrix (condensation and removal term)
%   B : Coagulation formation matrix (set input as empty if no
%   coagulation) CELL ARRAY
%   C : Coagulation removal matrix
%   n_init : Initial size distribution for the problem
%   t_interval : first and last time index in the vector [0,tmax]
%   h : time step for Crank-Nicolson
%   C_in : Border condition
%   dC_in : derivative of border condition
%   Source : source term (independent of size distribution function)
%
%   OUTPUTS
%   n_FEM : time evolution of size distribution
%   tt : time vector for the evolution
%
%% Intitializations for FEM and time evolution
tt = t_interval(1):h:t_interval(2);
tt = tt(:);
n_original = n_init;
n_FEM = zeros(length(n_original),length(tt));

% If boundary conditions exists
if ~isempty(C_in)
    
    % Applying boundary conditions
    M(1,:) = [];
    K(1,:) = [];
    M_b = M(:,1);
    K_b = K(:,1);
    M(:,1) = [];
    K(:,1) = [];
    K_init = K;
    n_init(1) = [];
    n_FEM(2:end,1) = n_init;
    
else
    
    % Without boundary condition
    n_FEM(:,1) = n_init;
    K_init = K;
    
end

% If Coagulation is present
if ~isempty(B) && ~isempty(C)
    
    for jj = 1:length(B)
        
        if ~isempty(C_in)
            
            % Applying boundary conditions for the FEM matrices
            BB = B{jj};
            CC = C{jj};
            BB(1,:) = [];
            BB(:,1) = [];
            CC(1,:) = [];
            CC(:,1) = [];
            
            % Formation of coagulation FEM matrix with boundary condition
            if jj > 1
                
                BC{jj-1,1} = BB-CC;
                
            end
            
        else
            
            % Formation of coagulation FEM matrix without boundary
            % conditions
            
            BC{jj,1} = B{jj}-C{jj};
            
        end
        
    end
    
else
    
    
    
    
end

% Calculating time evolutino for FEM
tic
for jj = 1:length(tt)-1
    
    % Check for if coagulation is present
    if exist ('BC','var')
        % Set K matrix to empty
        K = [];
        HH = [];
        
        % Linearizing cougulation FEM matrix with or without boundary
        % conditions
        if ~isempty(C_in)
            
            for ii = 1:length(BC)
                
                HH(ii,:) = n_FEM(2:end,jj)'*BC{ii};
                
            end
        else
            
            for ii = 1:length(BC)
                
                HH(ii,:) = n_FEM(:,jj)'*BC{ii};
                
            end
            
        end
        
        % Forming new stiffness matrix
        K = K_init+HH;
        
    end
    % Left Hand side for Crank-Nicolson method
    S_l = [eye(size(M))-(h/2).*(M\(K))];
    % Right hand side for Crank-Nicolson (This part alway exists)
    S_r = [eye(size(M))+(h/2).*(M\(K))];
    
    % If Source term is included
    if nargin > 9
        S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))))+source(tt(jj+1));
        S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))))+source(tt(jj));
        
        % Without source term
    elseif ~isempty(C_in) && nargin <= 9
        S2= (M\(-M_b*dC_in(tt(jj+1))+K_b*C_in(tt(jj+1))));
        S1= (M\(-M_b*dC_in(tt(jj))+K_b*C_in(tt(jj))));
        
        % Only source term but no boundary conditions
    elseif isempty(C_in) && nargin > 9
        S2= source(tt(jj+1));
        S1= source(tt(jj));
    end
    if exist('S1','var')
        n_FEM(2:end,jj+1) = S_l\S_r*n_FEM(2:end,jj) + (h/2)*S_l\(S1+S2);
        n_FEM(1,jj+1) = C_in(tt(jj+1));
%         n_FEM(n_FEM(:,jj+1)<0)=0;
    else
        n_FEM(:,jj+1) = S_l\(S_r)*n_FEM(:,jj);
%         n_FEM(n_FEM(:,jj+1)<0)=0;
    end
    
    % OPTIONAL
    % Positivity constraint
%     n_FEM((n_FEM(:,jj+1) < 0),jj+1) = 0;
    
    
    
end
evolution_time = toc;


end

