function [N_diff,n_diff,d_middle,d_edges,tt] = discrete_GDE_solver_ver2(d_min,d_max,coag_kernel,...
    n_init,source,time_step,t_max,solver,dN0,n0)
%DISCRETE_COAGULATION_SOLVER_VER2 calculates time evolution of coagulation
%equation using the discretized formulation by using correct formulation. This is valid only for
%evenly spaced bins. In addition, every bin have multiple of the monomer size.
% Time evolution method can be chosen between euler and crank-nicolson method.
% Note: Can be used either to particle volume v or particle diameter d.
%
% INPUTS:
%   d_min : middle point of smallest bin
%   d_max : maximum size for middle of the bin
%   coag_kernel : coagulation kernel function (must accept vectors as
%                 input)
%   cond_speed : condensation speed. In discrete formulation condensation
%                is assumed to be coagulation of smallest particles with bigger
%                particles. Theoretically number of particles in the firs
%                bin.
%   n_init : initial distribution for size distribution function
%            (N_i = n_i*Dd_i)
%   time_step : time step for the time evolution method
%   t_max : time for the time evolution
%   solver : decision to use either Euler or Crank-Nicolson method.
%   dN0 : is volume/diameter of monomer
%   n0 : Value for the size distribution function for monomer concentration
%
% OUTPUTS:
%   N_diff : time evolution of number distribution
%   n_diff : time evolution of size distribution
%   d_middle : middle points of bins
%   d_edges : eddges of the each bin
%   tt : time instances

% Making discretization for middle of bins
d_middle = [d_min:dN0:d_max]'; % Note: monomer is now used in discretization
d_width = d_middle(3)-d_middle(2); % Should be same for smallest size too?
d_edges = [d_middle-0.5*d_width,d_middle+0.5*d_width];
N0 = n0*d_width; % Value is changed for the number of particles

%% Coagulation formation and removal matrix creations
F = zeros(length(d_middle));
R = zeros(length(d_middle));

%These have to be done only ones
% index1 = find(d_middle == 2*d_middle(1));
[~,index1]=min(abs(d_middle-2*d_middle(1)));
index2 = find(d_middle == d_middle(end)-d_middle(1));

for ii = 1:index2
    
    F(index1+ii-1:end,ii) = coag_kernel(d_middle(ii),d_middle(1:index2+1-ii));
    
end
for ii = 1:length(d_middle)
    
    
    %     if d_middle(ii) >= 2*d_middle(1)
    %         F(ii:end,ii-1) = coag_kernel(d_middle(1:(end-ii+1),1),d_middle(ii-1));
    %     end
    
    R(:,ii) = coag_kernel(d_middle,d_middle(ii));
    
end

% Using the initial size distribution function to make the distribution
n_diff(:,1) = n_init(d_middle);
N_diff(:,1) = d_width*n_diff;

tt = 0:time_step:t_max;
% if (solver == 'Euler' | solver == 'euler')
if strcmp('Euler',solver) | strcmp('euler',solver)
    
    %     apu3 = repmat(d_middle(1:end-1,1),1,length(d_middle)-1);
    %     [apu3,I] = sort(tril(flipud(apu3)));
    for ii = 1:length(tt)-1
        %         apu1 = zeros(length(d_middle));
        %         apu4 = zeros(length(d_middle)-1);
        apu2 = repmat(N_diff(:,ii),1,length(N_diff(:,ii)));
        %         apu4 = tril(flipud(repmat(N_diff(1:end-1,ii),1,length(d_middle)-1)));
        %         apu4 = apu4(I);
        %         apu1(2:end,1:end-1) = apu4;
        
        apu1 = zeros(size(F));
        for jj = 1:length(N_diff(:,ii))-1
            
            apu1(jj+1:length(N_diff(:,ii)),jj) = N_diff(1:end-jj,ii);
            
        end
        
        % Condensation matrix (ditriagonal matrix)
        C = diag(-N0.*coag_kernel(dN0,d_middle))+diag(N0.*coag_kernel(dN0,d_middle(1:end-1)),-1);
        % Source term
        S = [source(tt(ii)),zeros(1,length(d_middle)-1)]';
        
        
        N_diff(:,ii+1) = N_diff(:,ii)+time_step*(0.5*F.*apu1-R.*apu2+C)*N_diff(:,ii)+time_step*S;
        n_diff(:,ii+1) = N_diff(:,ii+1)/d_width;
    end
    
    
elseif (strcmp('Crank-Nicolson',solver) | strcmp('crank-nicolson',solver) ...
        | strcmp('crank',solver) | strcmp('cranknicolson',solver))
    %     apu3 = repmat(d_middle(1:end-1,1),1,length(d_middle)-1);
    %     [apu3,I] = sort(tril(flipud(apu3)));
    
    for ii = 1:length(tt)-1
        %         apu1 = zeros(length(d_middle));
        %         apu4 = zeros(length(d_middle)-1);
        apu2 = repmat(N_diff(:,ii),1,length(N_diff(:,ii)));
        %         apu4 = tril(flipud(repmat(N_diff(1:end-1,ii),1,length(d_middle)-1)));
        %         apu4 = apu4(I);
        %         apu1(2:end,1:end-1) = apu4;
        apu1 = zeros(size(F));
%         for jj = 1:length(N_diff(:,ii))-index1+1
%             
%             apu1(jj+1:length(N_diff(:,ii)),jj) = N_diff(1:end-jj,ii);
%             
%         end
        for jj = 1:index2
            
            apu1(index1+jj-1:end,jj) = N_diff(1:index2+1-jj,ii);
            
        end
        
        % Condensation matrix (ditriagonal matrix)
        C = diag(-N0.*coag_kernel(dN0,d_middle))+diag(N0.*coag_kernel(dN0,d_middle(1:end-1)),-1);
        % Source term
        So = [source(tt(ii)),zeros(1,length(d_middle)-1)]';
        
        S = eye(length(d_middle))+(time_step/2)*(0.5*F.*apu1-R.*apu2+C);
        S2 = eye(length(d_middle))-(time_step/2)*(0.5*F.*apu1-R.*apu2+C);
        N_diff(:,ii+1) = (S2\S)*N_diff(:,ii)+S2\(time_step*So);
        n_diff(:,ii+1) = N_diff(:,ii+1)/d_width;
    end
    
    
else
    
    disp('ERROR. Unknown solver. Try Euler or Crank-Nicolson')
    
    
end





end

