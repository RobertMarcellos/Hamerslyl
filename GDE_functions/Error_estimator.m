function [relative_error,average_relative_error] = Error_estimator(disc_est,n_est,disc_acc,n_acc,disc_edges,err_mod)
%ERROR_ESTIMATOR calculates relative approximation error for the estimation
%methods by interpolating the approximation into the same size as the
%accurate solution.
%
% Inputs:
% disc_est  :   Elements and node points
% n_est     :   Approximation for the time evolution of the aerosol size
%               distribution
% disc_acc  :   discretization for the accurate solution (analytical or
%               discrete)
% n_acc     :   Analytical solution for GDE or the discrete GDE
% disc_edges:   Edges for the bins (for sectional method only)
% err_mod   :   Decision if the bin model error is calculated with linear of
%               histogramic interpolation
%
% Outputs:
% relative_error            : Relative error between the approximation and
%                             accurate solution
% average_relative_error    : Average from previous error
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2020
n_acc(isnan(n_acc)) = 0;

% Histogram error for the sectional method
if nargin == 6 && err_mod == 2
    
    n_disc_int = zeros(size(n_acc));
    
    % Outer iteration for time steps. Inner for the size axis
    for jj = 1:size(n_acc,2)
        for ii = 1:length(disc_est)
            index = [];
            index = find(disc_edges(ii) <= disc_acc & disc_edges(ii+1) >= disc_acc);
            
            n_disc_int(index,jj) = n_est(ii,jj);
            
        end
    end
    
    relative_error = 100*trapz(disc_acc,abs(n_disc_int-n_acc))./trapz(disc_acc,n_acc);
    average_relative_error = mean(relative_error);
    
elseif nargin == 6 && err_mod == 1
    
    % Linear error estimate for the sectional method
    temp1 = [];
    disc_est = disc_est(:);
    disc_est = [disc_edges(1);disc_est;disc_edges(end)];
    n_est = [n_est(1,:); n_est; n_est(end,:)];
    temp1 = abs(interp1(disc_est,n_est,disc_acc)-n_acc);
    temp1(isnan(temp1))=0;
    
    % Calculating error with the Eq. 24
    relative_error = 100*trapz(disc_acc,abs(temp1))./trapz(disc_acc,n_acc);
    average_relative_error = mean(relative_error);
    
    
else
    % Error for the FEM
    
    temp1 = [];
    temp1 = abs(interp1(disc_est,n_est,disc_acc)-n_acc);
    temp1(isnan(temp1))=0;
    
    relative_error = 100*trapz(disc_acc,abs(temp1))./trapz(disc_acc,n_acc);
    average_relative_error = mean(relative_error);
    
    
    
end

end

