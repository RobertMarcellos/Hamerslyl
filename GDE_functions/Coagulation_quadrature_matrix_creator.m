function [ B,C,B_petrov,C_petrov ] = Coagulation_quadrature_matrix_creator( nodes,kernel_function,gaussian_points,epsilon)
%COAGULATION_QUADRATURE_MATRIX_CREATOR calculates coagulation FEM matrixes
%by using the three or five point gaussian quadrature.

% INPUTS:
% nodes : node points for FEM solution
% kernel_function : coagulation kernel as a function (have to be apply to
%                   solve equation in a vector form)
% gaussian_points : Decision between 3 and 5 point Gaussian quadrature.
%                   Default is 3.
% epsilon : Upwind parameter for Petrov_Galerkin finite element method
%
% OUTPUTS:
% B : FEM matrix of the coagulation formation term
% C : FEM matrix of the coagulation loss term
% B_petrov : Petrov-Galerkin FEM matrix of the coagulation formation term
% C_petrov : Petrov-Galerkin FEM matrix of the coagulation loss term


nodes = nodes(:);
Vmin = min(nodes);
% Vmin = 0; % Testi

% Indexing different elements and getting value for the element boundaries
element_index = [1:length(nodes)-1;2:length(nodes)]';
element_boundary = [nodes(element_index(:,1)),nodes(element_index(:,2))];

% Element support as one basis function exist in two elements. Exception
% are the first and the last basis function which exist only in one
% element.
H_support = [nodes(1:end-2),nodes(2:end-1),nodes(3:end)];
H_support = [nan nodes(1) nodes(2);H_support; nodes(end-1),nodes(end),nan];

% Function works with only two inputs. Set the default to be 3 point
% gaussian quadrature.
if ~exist('gaussian_points','var') || isempty(gaussian_points)
    
    gaussian_points = 3;
    
end

% Initializing FEM matrixes
for l = 1:length(nodes)
    B{l,1} = sparse(zeros(length(nodes)));
    C{l,1} = sparse(zeros(length(nodes)));
    if nargin > 3
        B_petrov{l,1} = sparse(zeros(length(nodes)));
        C_petrov{l,1} = sparse(zeros(length(nodes)));
    end
end

% Looping through all the elements. Number of elements is the number of the
% nodes-1.
for jj = 1:length(nodes)-1 
    
    % Initializing local FEM formation matrix
    B_loc1 = zeros(length(nodes));
    B_loc2 = zeros(length(nodes)); 
    if nargin > 3
        
        B_PGloc1 = zeros(length(nodes));
        B_PGloc2 = zeros(length(nodes));
        
    end
    
    % Index for the current element
    e_ind = element_index(jj,:);
    
    % Looping all previous elements to end of the current element jj
    for kk = 1:min(e_ind(2),length(nodes)-1)
        
        
        % Finding the elements which produces particles corresbonding the
        % test function element phi_j. Always using the smaller element for
        % discretization so excluding the results which are smaller than
        % values in the element kk.
        [row,~] = find(H_support >= element_boundary(jj,1)-element_boundary(kk,2) ...
            & H_support <= element_boundary(jj,2)-element_boundary(kk,1));
        row = sort(unique(row));
        row(row < kk) = [];
        
        
        if ~isempty(row)
            
            % Looping through found basis functions phi_i
            for ii = 1:length(row)
                
                % Upper limit of the integration for phi_k can only be
                % v_j-Vmin or the v_k (end point of current element)
                up_lim = min(element_boundary(kk,2),element_boundary(jj,2)-Vmin); % Integration limit for the inner integration
                
                % Integration interval
                int_interval = [element_boundary(kk,1),up_lim];
              
                % Checking which boundaries of element jj exist in the
                % interval corresponding the summation of basis function
                % phi_i and the border of the element kk
                [~,col] = find(element_boundary(jj,:) >= H_support(row(ii),1)...
                    +int_interval(1) & element_boundary(jj,:) <= H_support(row(ii),3)+int_interval(2));
                interval_check = zeros(1,2);
                interval_check(col) = 1;
                
                
                % Different poosible integrations depending of the result
                % of the previous check. Setting limits for the outer
                % integral using these
                if (interval_check(1) ~= 0 || interval_check(2) ~= 0) || ...
                        (interval_check(1) == 0 && interval_check(2) == 0 && element_boundary(jj,2) >= int_interval(2)+H_support(row(ii),3) && element_boundary(jj,1) <= int_interval(1)+H_support(row(ii),1))
                    
                    % Both borders of element jj are between integration
                    % limits.
                    if interval_check == [1,1]
                        x1 = element_boundary(jj,1)-(H_support(row(ii),1)+element_boundary(kk,1));
                        x2 = (H_support(row(ii),3)+int_interval(2))-element_boundary(jj,2);
                        i_interval = [H_support(row(ii),1)+x1,H_support(row(ii),3)-x2];
                        v_interval = [element_boundary(jj,1),element_boundary(jj,2)];
                    
                        % Only the lower limit is
                    elseif interval_check == [1,0]
                        
                        x1 = element_boundary(jj,1)-(H_support(row(ii),1)+element_boundary(kk,1));
                        i_interval = [H_support(row(ii),1)+x1,H_support(row(ii),3)];
                        v_interval = [element_boundary(jj,1),H_support(row(ii),3)+int_interval(2)];
                        
                        % Only the upper limit is
                    elseif interval_check == [0,1]
                        
                        x2 = (H_support(row(ii),3)+int_interval(2))-element_boundary(jj,2);
                        i_interval = [H_support(row(ii),1),H_support(row(ii),3)-x2];
                        v_interval = [H_support(row(ii),1)+element_boundary(kk,1),element_boundary(jj,2)];
                    end
                    
                    
                    % Initializing gaussian quadrature for the outer
                    % integral
                    h_j = v_interval(2)-v_interval(1);
                    % Five point Gaussian quadrature
                    if gaussian_points == 5
                        w_j = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900]';
                        xi_j = [-(1/3)*sqrt(5+2*sqrt(10/7)),-(1/3)*sqrt(5-2*sqrt(10/7)),0,...
                            (1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))]';
                        
                        % Three point Gaussian quadrature
                    else
                        w_j = [5/9 8/9 5/9]';
                        xi_j = [-sqrt(3/5),0,sqrt(3/5)]';
                    end
                    
                    % Mapping to the master element [-1,1]
                    v_j = (0.5*h_j)*(xi_j+1) + v_interval(1);
                    phi_j1 = 0.5*(1-xi_j);
                    phi_j2 = 0.5*(1+xi_j);
                    
                    if nargin > 3
                        phiPG_j1 = 0.5*(1-xi_j)-(3/2)*epsilon*((4/h_j^2)*(v_j-v_interval(1)).*(-v_j+v_interval(2)));
                        phiPG_j2 = 0.5*(1+xi_j)+(3/2)*epsilon*((4/h_j^2)*(v_j-v_interval(1)).*(-v_j+v_interval(2)));
                        
                    end
                    
                    
                    % Integration corresponding the outer discretization
                    % using the gaussian quadrature function (at the end of
                    % the code)
                    for hh = 1:length(v_j)
                        % Check if the middle element of the basis function
                        % phi_i is in the integration interval. If it is,
                        % then the integration is done in two parts.
                        
                        if int_interval(1) <= v_j(hh)-H_support(row(ii),2) && int_interval(2) >= v_j(hh)-H_support(row(ii),2)
                            int_interval_1 = [int_interval(1), v_j(hh)-H_support(row(ii),2)];
                            int_interval_2 = [v_j(hh)-H_support(row(ii),2), int_interval(2)];
                            inner_integral1(hh,1) = quadrature(int_interval_1,v_j(hh),1,row(ii))+quadrature(int_interval_2,v_j(hh),1,row(ii));
                            inner_integral2(hh,1) = quadrature(int_interval_1,v_j(hh),2,row(ii))+quadrature(int_interval_2,v_j(hh),2,row(ii));
                        else
                            inner_integral1(hh,1) = quadrature(int_interval,v_j(hh),1,row(ii));
                            inner_integral2(hh,1) = quadrature(int_interval,v_j(hh),2,row(ii));
                        end 
                        
                    end
                    %
                    B_loc1(kk,row(ii)) = B_loc1(kk,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral1.*phi_j1);
                    %                 B_loc1(kk,row(ii)+1) = B_loc1(kk,row(ii)+1)+h_j*sum(w_j.*inner_integral12.*phi_j1);
                    B_loc1(kk+1,row(ii)) = B_loc1(kk+1,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral2.*phi_j1);
                    %                 B_loc1(kk+1,row(ii)+1) = B_loc1(kk+1,row(ii)+1)+h_j*sum(w_j.*inner_integral22.*phi_j1);
                    B_loc2(kk,row(ii)) = B_loc2(kk,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral1.*phi_j2);
                    %                 B_loc2(kk,row(ii)+1) = B_loc2(kk,row(ii)+1)+h_j*sum(w_j.*inner_integral12.*phi_j2);
                    B_loc2(kk+1,row(ii)) = B_loc2(kk+1,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral2.*phi_j2);
                    %                 B_loc2(kk+1,row(ii)+1) = B_loc2(kk+1,row(ii)+1)+h_j*sum(w_j.*inner_integral22.*phi_j2);
                    B_loc1(row(ii),kk) = B_loc1(kk,row(ii));
                    %                 B_loc1(row(ii)+1,kk) = B_loc1(kk,row(ii)+1);
                    B_loc1(row(ii),kk+1) = B_loc1(kk+1,row(ii));
                    %                 B_loc1(row(ii)+1,kk+1) = B_loc1(kk+1,row(ii)+1);
                    B_loc2(row(ii),kk) = B_loc2(kk,row(ii));
                    %                 B_loc2(row(ii)+1,kk) = B_loc2(kk,row(ii)+1);
                    B_loc2(row(ii),kk+1) = B_loc2(kk+1,row(ii));
                    %                 B_loc2(row(ii)+1,kk+1) = B_loc2(kk+1,row(ii)+1);
                    
                    if nargin > 3
                        
                        B_PGloc1(kk,row(ii)) = B_PGloc1(kk,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral1.*phiPG_j1);
                        B_PGloc1(kk+1,row(ii)) = B_PGloc1(kk+1,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral2.*phiPG_j1);
                        B_PGloc2(kk,row(ii)) = B_PGloc2(kk,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral1.*phiPG_j2);
                        B_PGloc2(kk+1,row(ii)) = B_PGloc2(kk+1,row(ii))+0.5*h_j/2*sum(w_j.*inner_integral2.*phiPG_j2);
                        B_PGloc1(row(ii),kk) = B_PGloc1(kk,row(ii));
                        B_PGloc1(row(ii),kk+1) = B_PGloc1(kk+1,row(ii));
                        B_PGloc2(row(ii),kk) = B_PGloc2(kk,row(ii));
                        B_PGloc2(row(ii),kk+1) = B_PGloc2(kk+1,row(ii));
                        
                    end
                end
            end
            
        end
        
    end
    
    % Placing the local FEM matrixes in to the 3D FEM matrix
    B{jj,1} = sparse(B{jj,1}+B_loc1);
    B{jj+1,1} = sparse(B{jj+1}+B_loc2);
    
    if nargin > 3
        B_petrov{jj,1} = sparse(B_petrov{jj,1}+B_PGloc1);
        B_petrov{jj+1,1} = sparse(B_petrov{jj+1}+B_PGloc2);
    end
    
    % Coagulation matrix C calculations
    % Handling always one element at the time so we need only to loop k
    % values

    % Initialization for Gaussian quadrature for the loss FEM matrix C
    
    % Width of the current element
    h_j = nodes(e_ind(2))-nodes(e_ind(1));
    % Five point Gaussian quadrature
    if gaussian_points == 5
        w_j = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900]';
        xi_j = [-(1/3)*sqrt(5+2*sqrt(10/7)),-(1/3)*sqrt(5-2*sqrt(10/7)),0,...
            (1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))]';
        % Three point Gaussian quadrature
    else
        w_j = [5/9 8/9 5/9]';
        xi_j = [-sqrt(3/5),0,sqrt(3/5)]';
    end
    
    
    v_j = 0.5*h_j*(xi_j+1) + nodes(e_ind(1));
    phi_j1 = 0.5*(1-xi_j);
    phi_j2 = 0.5*(1+xi_j);
    
    C_loc1 = zeros(length(nodes));
    C_loc2 = zeros(length(nodes));
    
    if nargin > 3
        phiPG_j1 = 0.5*(1-xi_j)-(3/2)*epsilon*((4/h_j^2)*(v_j-nodes(e_ind(1))).*(-v_j+nodes(e_ind(2))));
        phiPG_j2 = 0.5*(1+xi_j)+(3/2)*epsilon*((4/h_j^2)*(v_j-nodes(e_ind(1))).*(-v_j+nodes(e_ind(2))));
        
        C_PGloc1 = zeros(length(nodes));
        C_PGloc2 = zeros(length(nodes));
        
    end
    
    % Looping the inner integral values
    for kk = 1:length(nodes)-1
        
        % Initialization for Gaussian quadrature of the inner integral.
        ek_ind = element_index(kk,:);
        h_k = nodes(ek_ind(2))-nodes(ek_ind(1));
        v_kk = (h_k/2)*(xi_j+1) + nodes(ek_ind(1));
        phi_k1 = 0.5*(1-xi_j);
        phi_k2 = 0.5*(1+xi_j);
        
        for hh = 1:length(v_j)
            
            inner_integral1(hh,1) = h_k/2*sum(w_j.*kernel_function(v_j(hh),v_kk).*phi_k1);
            inner_integral2(hh,1) = h_k/2*sum(w_j.*kernel_function(v_j(hh),v_kk).*phi_k2);
            
        end
        
        % Forming the local FEM matrix
        C_loc1(jj,kk) = C_loc1(jj,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j1.*phi_j1);
        C_loc1(jj+1,kk) = C_loc1(jj+1,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j2.*phi_j1);
        C_loc1(jj,kk+1) = C_loc1(jj,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j1.*phi_j1);
        C_loc1(jj+1,kk+1) = C_loc1(jj+1,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j2.*phi_j1);
        
        
        C_loc2(jj,kk) = C_loc2(jj,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j1.*phi_j2);
        C_loc2(jj+1,kk) = C_loc2(jj+1,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j2.*phi_j2);
        C_loc2(jj,kk+1) = C_loc2(jj,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j1.*phi_j2);
        C_loc2(jj+1,kk+1) = C_loc2(jj+1,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j2.*phi_j2);
        
        if nargin > 3
            
            C_PGloc1(jj,kk) = C_PGloc1(jj,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j1.*phiPG_j1);
            C_PGloc1(jj+1,kk) = C_PGloc1(jj+1,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j2.*phiPG_j1);
            C_PGloc1(jj,kk+1) = C_PGloc1(jj,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j1.*phiPG_j1);
            C_PGloc1(jj+1,kk+1) = C_PGloc1(jj+1,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j2.*phiPG_j1);
            
            
            C_PGloc2(jj,kk) = C_PGloc2(jj,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j1.*phiPG_j2);
            C_PGloc2(jj+1,kk) = C_PGloc2(jj+1,kk)+h_j/2*sum(w_j.*inner_integral1.*phi_j2.*phiPG_j2);
            C_PGloc2(jj,kk+1) = C_PGloc2(jj,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j1.*phiPG_j2);
            C_PGloc2(jj+1,kk+1) = C_PGloc2(jj+1,kk+1)+h_j/2*sum(w_j.*inner_integral2.*phi_j2.*phiPG_j2);
            
        end

    end

    % Assigning local FEM matrix into 3D FEM matrix
    C{jj,1} = sparse(C{jj,1}+C_loc1);
    C{jj+1,1} = sparse(C{jj+1}+C_loc2);
    
    if nargin > 3
        
        C_petrov{jj,1} = sparse(C_petrov{jj,1}+C_PGloc1);
        C_petrov{jj+1,1} = sparse(C_petrov{jj+1}+C_PGloc2);
        
    end
    
    
    
    
    
end

% Gaussian quadrature function for the inner integrals for formation term.
    function int_value = quadrature(integration_limits,v_j_value,phi_k,basis_i)
        % integration_limits : Limits for the current integration
        % v_j_value : current value for outer integral discretization
        % phi_k : decision if the line is inqreasing or decreasing in the
        %         current element
        % basis_i : Current basis function phi_i
        
        % Refining the integration limits so that all the basis functions
        % get values in the integration.
        integration_limits(2) = min([integration_limits(2),v_j_value-min(H_support(row(ii),:)),v_j_value-Vmin]);
        integration_limits(1) = max([integration_limits(1),v_j_value-max(H_support(row(ii),:))]);
        
        
        
        % Initializing gaussian quadrature
        h = integration_limits(2)-integration_limits(1);
        % Five point gaussian quadrature
        if gaussian_points == 5
                    w = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900]';
                    xi = [-(1/3)*sqrt(5+2*sqrt(10/7)),-(1/3)*sqrt(5-2*sqrt(10/7)),0,...
                        (1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))]';
            % Three point Gaussian quadrature
        else
            w = [5/9 8/9 5/9]';
            xi = [-sqrt(3/5),0,sqrt(3/5)]';
        end
        v_k = 0.5*h*(xi+1) + integration_limits(1);
        if phi_k == 1
            phi = 0.5*(1+xi);
        elseif phi_k == 2
            phi = 0.5*(1-xi);
        end
        
        int_value = (h/2)*sum(w.*kernel_function(v_j_value-v_k,v_k).*phi...
            .*basis_fun(nodes,basis_i,v_j_value-v_k,min(nodes),max(nodes)));
        
    end

% Finite element basis function. Generally, these exists in two element and
% this function is used to get corresponding value for the gaussian
% quadratures.
    function f = basis_fun(node_points,num_of_node,x,Vmin,limit)
        f = zeros(size(x));
        g = node_points;
        k = num_of_node;
        % x = unique(x);
        
        % First and last element have to be calculated separately because there
        % are not whole triangles.
        if k == 1
            %     interval = [-g(n+1) g(n+1)];
            interval = [g(k) g(k+1)];
            % Constant for increasing and decreasing line
            b2 = 1+1/(interval(2)-g(k))*g(k);
            % Checking which points are at the increasing line
            ind = find(x < g(k) & x >= interval(1));
            %     f(ind) = (1/(g(k)-interval(1)).*x(ind)+b1);
            f(ind) = 0.*x(ind); % Setting these to be zero because no increasing line
            
            % f(x > g(k) & x <= interval(2)) = (-1/(interval(2)-g(k))*(x)+b2);
            % Indexes in decreasing line
            ind2 = find(x >= g(k) & x <= interval(2));
            f(ind2) = (-1/(interval(2)-g(k)).*x(ind2)+b2);
            
            % Last node accordingly to the first node
        elseif k == length(g)
            
            %     interval = [g(k-1) g(k)+abs(g(k)-g(k-1))];
            interval = [g(k-1) g(k)];
            b1 = 1-1/(g(k)-interval(1))*g(k);
            b2 = 1+1/(interval(2)-g(k))*g(k);
            ind = find(x <= g(k) & x >= interval(1));
            f(ind) = (1/(g(k)-interval(1)).*x(ind)+b1);
            % f(x > g(k) & x <= interval(2)) = (-1/(interval(2)-g(k))*(x)+b2);
            ind2 = find(x > g(k) & x <= interval(2));
            f(ind2) = (-1/(interval(2)-g(k)).*x(ind2)+b2);
            
        else
            interval = [g(k-1) g(k+1)];
            b1 = 1-1/(g(k)-interval(1))*g(k);
            b2 = 1+1/(interval(2)-g(k))*g(k);
            ind = find(x <= g(k) & x >= interval(1));
            f(ind) = (1/(g(k)-interval(1)).*x(ind)+b1);
            % f(x > g(k) & x <= interval(2)) = (-1/(interval(2)-g(k))*(x)+b2);
            ind2 = find(x > g(k) & x <= interval(2));
            f(ind2) = (-1/(interval(2)-g(k)).*x(ind2)+b2);
        end
        % f(isnan(f)) = 1;
        
        % If there are 5 inputs this if conditions sets the larger values than the
        % limit value to be zero
        if nargin == 5
            f(x > limit) = 0;
        end
        % Also for the lower end
        f(x < Vmin) = 0;
    end
% Petrov-Galerkin test function

end

