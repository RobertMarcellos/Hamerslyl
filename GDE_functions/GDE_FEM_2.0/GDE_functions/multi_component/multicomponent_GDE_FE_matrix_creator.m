function [A,G1,G2,B,C,R,A_PG,G1_PG,G2_PG,B_PG,C_PG,R_PG] = multicomponent_GDE_FE_matrix_creator(nodes,growth_rate,dgrowth_rate,coag_kernel,epsilon,removal)
%MULTICOMPONENT_GDE_FE_MATRIX_CREATOR creates FE matrices with standard test function and with Petrov-Galerkin test function. 
% This function can be used to form needed matrices
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2022
%
% Inputs:
% nodes : node points for the FE elements
% growth_rate : cell array containing growth rates for the different
% compounds
% dgrowth_rate : cell array containing the derivatives of the growth rates
%                if set empty, boundary condition should be applied in time
%                evolution calculations
%                NOTE: It is preferable to set boundary condition
% coag_kernel : coagulation kernel function
% epsilon : Petrov-Galerkin parameter for the PGFEM
% removal : function for the removal of the particles (must be time independent)
%
% Outputs:
% A : mass matrix for FEM
% G1 : cell array containing parts of the stiffness matrix
% G2 : more parts of the stiffness matrices
% B : matrix for the coagulation formation term
% C : matrix for the coagulation loss term
% R : Matrix for the removal term
% REST parameters are for PGFEM

% Points for the gaussian quadrature
ksi = [1/2-sqrt(3)/6 , 1/2+sqrt(3)/6]';
w = [1/2 1/2]';


% Basis functions
phi1 = 1-ksi;
phi2 = ksi;
dphi1 = [-1 -1]';
dphi2 = [1 1]';

g = nodes(:);
H = [1:length(g)-1;2:length(g)]';


% Initializations for the FE matrices
[rg,sg]=size(g); 
[rH,sH]=size(H);
G2 = zeros(rg,rg); 
A = zeros(rg,rg);
G2_PG = zeros(rg,rg);
A_PG = zeros(rg,rg);
R = zeros(rg,rg);
R_PG = zeros(rg,rg);

% Renaming the functions if included
if nargin > 5 && ~isempty(removal)
    rfun = removal;
end

% Setting epsilon to empty array if the number of inputs is less than 5
if nargin < 5

    epsilon = [];

end

% Petrov-Galerkin basis functions
if  nargin > 4 && ~isempty(epsilon)
    PGphi1 = epsilon*(6*ksi.^2-6*ksi);
    PGphi2 = epsilon*(-6*ksi.^2+6*ksi);

    d_PGphi1 = epsilon*(12*ksi-6);
    d_PGphi2 = epsilon*(-12*ksi+6);
end



% Looping through elements of the FE matrices
for ii=1:rH
    ind=H(ii,:); 
    gind=g(ind,:);

    h = gind(2)-gind(1);
    x = h*ksi;

    % Massamatriisi
    m = [h*sum(w.*phi1.*phi1) h*sum(w.*phi1.*phi2);
        h*sum(w.*phi2.*phi1) h*sum(w.*phi2.*phi2)];

    if nargin > 5 && ~isempty(removal)
        % Removal matrix (IF constant in time)
        rem = -[h*sum(w.*rfun(gind(1)+x).*phi1.*phi1) h*sum(w.*rfun(gind(1)+x).*phi1.*phi2);
            h*sum(w.*rfun(gind(1)+x).*phi2.*phi1) h*sum(w.*rfun(gind(1)+x).*phi2.*phi2)];
        R(ind,ind)=R(ind,ind)+(rem)';

    end

    % PGFEM matrix for the removal term
    if nargin > 5 && ~isempty(epsilon) &&  ~isempty(removal)

        rem_PG = -[h*sum(w.*rfun(gind(1)+x).*phi1.*PGphi1) h*sum(w.*rfun(gind(1)+x).*phi1.*PGphi2);
            h*sum(w.*rfun(gind(1)+x).*phi2.*PGphi1) h*sum(w.*rfun(gind(1)+x).*phi2.*PGphi2)];
        R_PG(ind,ind)=R_PG(ind,ind)+(rem+rem_PG)';

    end

    % Sum of growth rate functions
    H_sum = zeros(size(x));
    dH_sum = zeros(size(x));

    % Forming the needed matrices for the gaussian quadrature
    for jj = 1:length(growth_rate)
        fun = growth_rate{jj};
        if ~isempty(dgrowth_rate)

            dfun = dgrowth_rate{jj};
            dH_sum = dH_sum + dfun(gind(1)+x);

        end

        H_sum = H_sum + fun(gind(1)+x);



    end

    if ~isempty(dgrowth_rate)

        % Stiffness matrix for the FEM (Without border condition)
        k1 = -[h*sum(w.*phi1.*phi1.*H_sum) h*sum(w.*phi1.*phi2.*H_sum);
            h*sum(w.*phi2.*phi1.*H_sum) h*sum(w.*phi2.*phi2.*H_sum)];

        k2 = -[sum(w.*dphi1.*phi1.*(gind(1)+x).*H_sum) sum(w.*dphi1.*phi2.*(gind(1)+x).*H_sum);
            sum(w.*dphi2.*phi1.*(gind(1)+x).*H_sum) sum(w.*dphi2.*phi2.*(gind(1)+x).*H_sum)];

        k3 = -[h.*sum(w.*phi1.*phi1.*(gind(1)+x).*dH_sum) h.*sum(w.*phi1.*phi2.*(gind(1)+x).*dH_sum);
            h.*sum(w.*phi2.*phi1.*(gind(1)+x).*dH_sum) h.*sum(w.*phi2.*phi2.*(gind(1)+x).*dH_sum)];

        % Assembling the local matrix
        k = k1 + k2 + k3;

    else % With the boundary condition

        k2 = [sum(w.*phi1.*dphi1.*(gind(1)+x).*H_sum) sum(w.*phi1.*dphi2.*(gind(1)+x).*H_sum);
            sum(w.*phi2.*dphi1.*(gind(1)+x).*H_sum) sum(w.*phi2.*dphi2.*(gind(1)+x).*H_sum)];

        k = k2;


    end



    % Inserting the local matrix into the global matrix
    A(ind,ind)=A(ind,ind)+(m)';
    G2(ind,ind) = G2(ind,ind)+k';

    if ~isempty(epsilon) && nargin > 4

        % Mass matrix for the PGFEM
        m_PG = [h*sum(w.*phi1.*PGphi1) h*sum(w.*phi1.*PGphi2);
            h*sum(w.*phi2.*PGphi1) h*sum(w.*phi2.*PGphi2)];

        if ~isempty(dgrowth_rate)

            % Stiffness matrix for the PGFEM (Without border condition)

            % Calculating the local matrix for the PGFEM
            k1_PG = -[h*sum(w.*phi1.*PGphi1.*H_sum) h*sum(w.*phi1.*PGphi2.*H_sum);
                h*sum(w.*phi2.*PGphi1.*H_sum) h*sum(w.*phi2.*PGphi2.*H_sum)];

            k2_PG = -[sum(w.*dphi1.*PGphi1.*(gind(1)+x).*H_sum) sum(w.*dphi1.*PGphi2.*(gind(1)+x).*H_sum);
                sum(w.*dphi2.*PGphi1.*(gind(1)+x).*H_sum) sum(w.*dphi2.*PGphi2.*(gind(1)+x).*H_sum)];

            k3_PG = -[h.*sum(w.*phi1.*PGphi1.*(gind(1)+x).*dH_sum) h.*sum(w.*phi1.*PGphi2.*(gind(1)+x).*dH_sum);
                h.*sum(w.*phi2.*PGphi1.*(gind(1)+x).*dH_sum) h.*sum(w.*phi2.*PGphi2.*(gind(1)+x).*dH_sum)];

            k_PG = k1 + k2 + k3 + k1_PG + k2_PG + k3_PG;

        else % Stiffness matrix for the PGFEM (Without boundary condition)


            k2_PG = [sum(w.*phi1.*d_PGphi1.*(gind(1)+x).*H_sum) sum(w.*phi1.*d_PGphi2.*(gind(1)+x).*H_sum);
                sum(w.*phi2.*d_PGphi1.*(gind(1)+x).*H_sum) sum(w.*phi2.*d_PGphi2.*(gind(1)+x).*H_sum)];
            %
            k_PG = k2 + k2_PG;

        end

        % Inserting the local matrix into the global matrix
        A_PG(ind,ind)=A_PG(ind,ind)+(m+m_PG)';
        G2_PG(ind,ind) = G2_PG(ind,ind)+k_PG';



    end


end

% Calculating the other stiffness matrix
G1 = cell(length(growth_rate),1);
G1_PG = cell(length(growth_rate),1);

for jj = 1:length(growth_rate)
    fun = growth_rate{jj};
    K1 = zeros(rg,rg);
    K1_PG = zeros(rg,rg);
    for ii=1:rH
        ind=H(ii,:); 
        gind=g(ind,:); 

        h = gind(2)-gind(1);
        x = h*ksi;

        k = [h*sum(w.*phi1.*phi1.*fun(gind(1)+x)) h*sum(w.*phi1.*phi2.*fun(gind(1)+x));
            h*sum(w.*phi2.*phi1.*fun(gind(1)+x)) h*sum(w.*phi2.*phi2.*fun(gind(1)+x))];
        K1(ind,ind)=K1(ind,ind)+k';

        % Petrov-Galerkin matrix
        if ~isempty(epsilon)
            ind=H(ii,:); 
            gind=g(ind,:); 

            h = gind(2)-gind(1);
            x = h*ksi;

            k_PG = [h*sum(w.*phi1.*PGphi1.*fun(gind(1)+x)) h*sum(w.*phi1.*PGphi2.*fun(gind(1)+x));
                h*sum(w.*phi2.*PGphi1.*fun(gind(1)+x)) h*sum(w.*phi2.*PGphi2.*fun(gind(1)+x))];

            % Assembling the global matrix
            K1_PG(ind,ind)=K1_PG(ind,ind)+(k+k_PG)';
        else

            K1_PG = [];
        end

    end
    G1{jj,1} = K1;
    G1_PG{jj,1} = K1_PG;

end

% Calculating coagulation matrices if the coagulation kernel is given.

% The case where also the PG parameter is given
if nargin > 3 && ~isempty(coag_kernel) && ~isempty(epsilon)

    % Coagulation FE matrices with 5-point quadrature
    [ B,C,B_PG,C_PG ] = Coagulation_MC_quadrature_matrix_creator(g,coag_kernel,5,epsilon);

    % Case where the parameter is not given
elseif nargin > 3 && ~isempty(coag_kernel) && isempty(epsilon)
    % Coagulation FE matrices with 5-point quadrature
    [ B,C] = Coagulation_MC_quadrature_matrix_creator(g,coag_kernel,5);
    B_PG = [];
    C_PG = [];
else
    B = [];
    B_PG = [];

    C = [];
    C_PG = [];

end


end
