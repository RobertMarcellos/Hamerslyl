function [m_diff,d_widths,C1,C2,R,X,Beta] = multicomponent_GDE_sectional_matrices_volume(interval,num_of_bins,growth_rate,removal,coag_kernel)
%MULTICOMPONENT_GDE_SECTIONAL_MATRICES calculates matrices for the
%sectional method. Condensation matrices are calculated with the 1st order
%difference. Coafulation matrices are calculated with the size splitting
%operator.
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2022
%
%   Explanation for inputs:
%   interval = [min,max] exponents for the size interval
%   num_of_bins = number of bins
%   growth_rate = cell array of the growth rates of the different
%                 substances
%   removal = removal rate of different substance
%   coag_kernel = coagulation kernel function
%
%   Explanation for outputs:
%
d_edges = logspace(interval(1),interval(2),num_of_bins)';
d_widths = diff(d_edges);
m_diff = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins

if ~isempty(growth_rate)
    
    C1 = cell(size(growth_rate));
    %     C2 = zeros(length(m_diff));
    
    H_tot = zeros(length(m_diff),1);
    for ii = 1:length(growth_rate)
        
        fun = growth_rate{ii};
        C1{ii} = diag(fun(m_diff));
        H_tot = H_tot + fun(m_diff);
        
        
    end
    
    bb = 1./diff(m_diff);
    bb = [1/m_diff(1); bb];
    G31 = -bb.*m_diff.*H_tot;
    G32 = [1;d_widths(2:end-1)./d_widths(1:end-2)].*bb(2:end).*m_diff(1:end-1).*H_tot(1:end-1);
    for ii = 1:size(G31,1)
        
        if ii == 1
            G_tot(ii,ii) = G31(ii);
        else
            
            G_tot(ii,ii-1:ii) = [G32(ii-1), G31(ii)];
            
        end
        
    end
    
    C2 = G_tot;
    
    
    
else
    
    C1 = [];
    C2 = [];
    
    
    
end

if  nargin > 3
    
    if ~isempty(removal)

        % If removal is different for each species
%         R = cell(size(removal));
%         
%         for ii = 1:length(removal)
%             
%             rfun = removal{ii};
%             R{ii} = diag(rfun(m_diff));
%             
%             
%         end
        R = diag(removal(m_diff));
        
        
    else
        
        R = [];
        
        
    end
    
    
    
else
    
    R = [];
    
    
    
end

if nargin > 4
    
    if ~isempty(coag_kernel)
        
        X = Size_splitting_operator(m_diff);
        for j = 1:length(m_diff)
            
            Beta(:,j) = coag_kernel(m_diff,m_diff(j));
            
        end
        
        
    else
        
        X = [];
        Beta = [];
        
        
    end
    
    
else
    
    X = [];
    Beta = [];
    
    
end

end

function [split_operator] = Size_splitting_operator(v_middle)
%SIZE_SPLITTING_OPERATOR function returns cell array containing  matrixes
%for particle formation in coagulation

%   Inputs:
%   v_middle = center volumes of different size bins
%
%   Outputs:
%   split_operator = Cell array containing size-splitting matrixes for different
%   size bins
split_operator = cell(length(v_middle),1);
for kk = 1:length(v_middle)
%     X_apu = [];
    X_apu = zeros(length(v_middle));
    
    for ii = 1:length(v_middle)
        
        for jj = 1:length(v_middle)
            
            if kk == 1
                if v_middle(ii)+v_middle(jj) >= v_middle(kk) && v_middle(ii)+v_middle(jj) < v_middle(kk+1)
                    
                    X_apu(ii,jj) = (v_middle(kk+1)-(v_middle(ii)+v_middle(jj)))/(v_middle(kk+1)-v_middle(kk));
                    
                else
                    
                    X_apu(ii,jj) = 0;
                    
                end
                
                
            elseif kk == length(v_middle)
%                 v_inter = 2*v_middle(kk)-v_middle(kk-1); % Interpolating extra point into the end to avoid excess loss
%                 
%                 if v_middle(ii)+v_middle(jj) >= v_middle(kk) && v_middle(ii)+v_middle(jj) < v_inter
%                     
%                     X_apu(ii,jj) = (v_inter-(v_middle(ii)+v_middle(jj)))/(v_inter-v_middle(kk));
                    
                if v_middle(ii)+v_middle(jj) >= v_middle(kk-1) && v_middle(ii)+v_middle(jj) < v_middle(kk)
                    
                    X_apu(ii,jj) = ((v_middle(ii)+v_middle(jj))-v_middle(kk-1))/(v_middle(kk)-v_middle(kk-1));
                    
                    
                else
                    
                    X_apu(ii,jj) = 0;
                    
                end
  
            else
                
                if v_middle(ii)+v_middle(jj) >= v_middle(kk) && v_middle(ii)+v_middle(jj) < v_middle(kk+1)
                    
                    X_apu(ii,jj) = (v_middle(kk+1)-(v_middle(ii)+v_middle(jj)))/(v_middle(kk+1)-v_middle(kk));
                    
                elseif v_middle(ii)+v_middle(jj) >= v_middle(kk-1) && v_middle(ii)+v_middle(jj) < v_middle(kk)
                    
                    X_apu(ii,jj) = ((v_middle(ii)+v_middle(jj))-v_middle(kk-1))/(v_middle(kk)-v_middle(kk-1));
                    
                    
                else
                    
                    X_apu(ii,jj) = 0;
                    
                end
                
                
            end
            
        end
    end
    
    split_operator{kk,1} = sparse(X_apu);
end
end

