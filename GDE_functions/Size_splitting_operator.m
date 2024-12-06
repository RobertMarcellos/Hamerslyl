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

% Size split operator corresponding the Eq. 20
for kk = 1:length(v_middle)
%     X_apu = [];
    X_apu = zeros(length(v_middle));
    
    for ii = 1:length(v_middle)
        
        for jj = 1:length(v_middle)
            
            % First and last bin are handled seperately
            if kk == 1
                if v_middle(ii)+v_middle(jj) >= v_middle(kk) && v_middle(ii)+v_middle(jj) < v_middle(kk+1)
                    
                    X_apu(ii,jj) = (v_middle(kk+1)-(v_middle(ii)+v_middle(jj)))/(v_middle(kk+1)-v_middle(kk));
                    
                else
                    
                    X_apu(ii,jj) = 0;
                    
                end
                
                
            elseif kk == length(v_middle)
                    
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