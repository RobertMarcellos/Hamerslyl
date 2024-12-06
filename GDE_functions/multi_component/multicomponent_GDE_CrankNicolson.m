function [q_tot,q_i,nodes,tt] = multicomponent_GDE_CrankNicolson(nodes,initial_dist,time_step,t_max,LHS,RHS1,RHS2,coag_form,coag_rem,removal_term,BC,dBC)
% multicomponent_GDE_CrankNicolson calculates the time evolution for volume
% concentrations q_i with Crank-Nicolson method (C-N method)
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2022
%
% INPUTS:
% Nodes         :   Node points distributed to chosen spatial interval
% initial_dist  :   Initial distributions of q_i as a cell array containing
%                   vector of initial distribution
% time_step     :   Time step for the C-N iteration
% t_max         :   Time instance for the end of iterations
% LHS           :   Finite element matrix of the left hand side of FEM. The
%                   dimension must be length(Nodes)*number of compounds
% RHS1          :   Cell array containing matrices which are individual for
%                   each component (multicomponent_GDE_FE_matrix_creator
%                   gives this as a output).
% RHS2          :   The FE matrix which is same for each compound
% coag_form     :   Coagulation formation FE matrix (set empty if no
%                   coagulation)
% coag_rem      :   Coagulation removal FE matrix (set empty if no
%                   coagulation)
% removal_term  :   Removal process as a function handle or as a individual
%                   matrix
% BC            :   Boundary condition for the FEM as a function
% dBC           :   Derivative (true or approximative) of previous function
%
% OUTPUTS:
% q_tot         :   Summation from volume concentrations of q_i
% q_i           :   Individual time evolutions for volume concentrations
%                   q_i
% nodes         :   Node points of FEM approximation
% tt            :   Vector containing time instances

% Matrix assembly
A = zeros(length(RHS1)*length(nodes));
G1 = zeros(length(RHS1)*length(nodes));
G2 = zeros(length(RHS1)*length(nodes));
Rem = zeros(length(RHS1)*length(nodes));

% Coagulation matrix assembly (if coagulation exists)
if ~isempty(coag_form)
    for ii = 1:length(coag_form)

        B_temp = coag_form{ii};
        C_temp = coag_rem{ii};

        BC_coag{ii,1} = sparse(B_temp-C_temp);


    end
end

% Creating other FEM matrices
for ii = 1:length(RHS1)

    A((ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes),(ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes)) = LHS;
    G2((ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes),(ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes)) = RHS2;

    % Applying removal process if one exists
    if nargin > 9
        if ~isempty(removal_term)

            if isa(removal_term,'function_handle')

                Rem((ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes),(ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes)) ...
                    = diag(removal_term(nodes));
            else

                Rem((ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes),(ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes)) ...
                    = removal_term;

            end

        end
    end
    for jj = 1:length(RHS1)

        G1((ii-1)*length(nodes)+1:(ii-1)*length(nodes)+length(nodes),(jj-1)*length(nodes)+1:(jj-1)*length(nodes)+length(nodes)) = RHS1{ii};

    end

end

if nargin > 9
    if ~isempty(removal_term)

        G2 = G2 - Rem;

    end
end

A = sparse(A);
G1 = sparse(G1);
G2 = sparse(G2);

% Applying possible boundary conditions and forming initial distributions
if nargin > 10

    % Collecting matrices for the boundary conditions
    M_b = length(RHS1)*zeros(length(nodes)-1,length(RHS1));
    K_b = length(RHS1)*zeros(length(nodes)-1,length(RHS1));
    K_b2 = length(RHS1)*zeros(length(nodes)-1,length(RHS1));

    for ii = 1:length(RHS1)


        M_b((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),ii) = LHS(2:end,1);
        if ~isempty(removal_term)
            if isa(removal_term,'function_handle')
                R_temp = diag(removal_term(nodes));
                K_b((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),ii) = RHS2(2:end,1)-R_temp(2:end,1);
            else
                K_b((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),ii) = RHS2(2:end,1)-removal_term(2:end,1);
            end

        else
            K_b((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),ii) = RHS2(2:end,1);
        end
        temp = RHS1{ii};
        for jj = 1:length(RHS1)
            K_b2((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),jj) = temp(2:end,1);
        end

    end

    q_tot_temp = [];
    for ii = 1:length(RHS1)

        temp = initial_dist{ii};
        temp(1) = [];
        q_tot_temp = [q_tot_temp;temp];

    end

else

    q_tot = [];
    for ii = 1:length(RHS1)

        temp = [];
        temp = initial_dist{ii};
        q_tot = [q_tot;temp];

    end

end

% Support for coagulation will be added later

tt = [0:time_step:t_max]';

% Removing additional lines if boundary condition is applied

if isempty(coag_form) && nargin > 10
    
    A_temp = A;
    G1_temp = G1;
    G2_temp = G2;
    
    A = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));

    G1 = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));
    G2 = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));
    
    for ii = 1:length(RHS1)
        for jj = 1:length(RHS1)
            
            A((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                A_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            G1((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                G1_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            G2((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                G2_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            
        end
    end
    

    S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2))];
    S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2))];
    inv_mat = S1\S2;

elseif nargin > 10
    
    A_temp = A;
    G1_temp = G1;
    G2_temp = G2;
    
    A = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));

    G1 = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));
    G2 = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));
    
    for ii = 1:length(RHS1)
        for jj = 1:length(RHS1)
            
            A((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                A_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            G1((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                G1_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            G2((ii-1)*(length(nodes)-1)+1:(ii)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                G2_temp((ii-1)*length(nodes)+2:(ii)*length(nodes),(jj-1)*length(nodes)+2:(jj)*length(nodes));
            
        end
    end



end

for ii = 1:length(tt)-1

    % Applying boundary condition
    if nargin > 10

        % If coagulation is accounted for
        if ~isempty(coag_form)

            HH = zeros(length(nodes));
            reuna = [];
            for kk = 1:length(BC)

                fun = BC{kk};

                reuna(:,kk) = [fun(ii);q_tot_temp((kk-1)*(length(nodes)-1)+1:(kk)*(length(nodes)-1),end)];

            end
            for jj = 1:length(nodes)

                HH(jj,:) = (BC_coag{jj,1}*(sum(reuna,2)))';

            end
            HH(1,:) = [];
            HH_BC = HH(:,1);
            HH(:,1) = [];
            temp = length(RHS1)*zeros(length(nodes)-1,length(RHS1));
            HH_block = zeros(length(RHS1)*length(nodes)-length(RHS1),length(RHS1)*length(nodes)-length(RHS1));
            for jj = 1:length(RHS1)
                temp((jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1),jj) = HH_BC;
                HH_block((jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1),(jj-1)*(length(nodes)-1)+1:(jj)*(length(nodes)-1)) = ...
                HH;
            end

            % Applying the effect of coagulation and taking account the
            % border condition
            K_b_new = K_b + temp;

           
            S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2+HH_block))];
            S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2+HH_block))];
            inv_mat = S1\S2;

            q_in_cur = zeros(length(BC),1);
            q_in_next = zeros(length(BC),1);
            dq_in_cur = zeros(length(BC),1);
            dq_in_next = zeros(length(BC),1);
            for kk = 1:length(BC)

                fun = BC{kk};
                dfun = dBC{kk};

                q_in_cur(kk,1) = fun(ii);
                q_in_next(kk,1) = fun(ii+1);

                dq_in_cur(kk,1) = dfun(ii);
                dq_in_next(kk,1) = dfun(ii+1);


            end


            q_tot_temp(:,ii+1) = inv_mat*q_tot_temp(:,ii) + 0.5*time_step*(S1\(A\K_b2))*(q_in_cur+q_in_next) + 0.5*time_step*(S1\(A\K_b_new))*(q_in_cur+q_in_next) - 0.5*time_step*(S1\(A\M_b))*(dq_in_cur+dq_in_next);

            % Positivity constraint
            q_tot_temp(q_tot_temp < 0) = 0;



            % Case if coagulation is not happening
        else

            q_in_cur = zeros(length(BC),1);
            q_in_next = zeros(length(BC),1);
            dq_in_cur = zeros(length(BC),1);
            dq_in_next = zeros(length(BC),1);
            for kk = 1:length(BC)

                fun = BC{kk};
                dfun = dBC{kk};

                q_in_cur(kk,1) = fun(ii);
                q_in_next(kk,1) = fun(ii+1);

                dq_in_cur(kk,1) = dfun(ii);
                dq_in_next(kk,1) = dfun(ii+1);


            end

            q_tot_temp(:,ii+1) = inv_mat*q_tot_temp(:,ii) + (time_step/2)*S1\(A\(K_b*(q_in_cur+q_in_next)+K_b2*(q_in_cur+q_in_next)-M_b*(dq_in_cur+dq_in_next)));


            % Positivity constraint
            q_tot_temp(q_tot_temp < 0) = 0;


        end

    else % No boundary condition (FEM should always have boundary condition)

        % With coagulation
        if ~isempty(coag_form)

            HH = zeros(length(q_tot(:,1))/2);
            HH_block = zeros(size(A));
            for jj = 1:length(nodes)

                HH(jj,:) = (BC_coag{jj,1}*(q_tot(1:size(q_tot,1)/2,end) + q_tot(size(q_tot,1)/2+1:end,end)))';

            end

            HH_block = blkdiag(HH,HH);
             if exist('removal_term','var')
                if ~isempty(removal_term) 
                    S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2+HH_block-Rem))];
                    S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2+HH_block-Rem))];
                else
                    S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2+HH_block))];
                    S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2+HH_block))];
                
                end
             else

                 S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2+HH_block))];
                 S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2+HH_block))];
                
             end
            inv_mat = S1\S2;
           
            q_tot(:,ii+1) = q_tot(:,ii+1)+time_step*(A\(G1+G2+HH_block)*q_tot(:,end));
            % Positivity constraint
            q_tot(q_tot < 0) = 0;


        else % Without coagulation

            if exist('removal_term','var')
                if ~isempty(removal_term) 
                    S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2-Rem))];
                    S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2-Rem))];
                else
                    S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2))];
                    S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2))];
                
                end
             else

                 S1 = [sparse(eye(size(A)))-(time_step/2)*(A\(G1+G2))];
                 S2 = [sparse(eye(size(A)))+(time_step/2)*(A\(G1+G2))];
                
             end
            
            inv_mat = S1\S2;

            q_tot(:,ii+1) = inv_mat*q_tot(:,ii);

            % Positivity constraint
            q_tot(q_tot < 0) = 0;

        end



    end


end

% Calculating different variables
if nargin > 10 % With boundary condition

    temp = size(q_tot_temp,1)/length(RHS1);

    for ii = 1:length(RHS1)

        % Applying boundary condition
        if ii == 1

            fun = BC{ii};
            source = zeros(1,length(tt));
            for jj = 1:length(tt)

                source(1,jj) = fun(jj);

            end

            %             source = BC{ii};
            q_i{ii} = [source ; q_tot_temp(1:length(nodes)-1,:)];

        else


            fun = BC{ii};
            source = zeros(1,length(tt));
            for jj = 1:length(tt)

                source(1,jj) = fun(jj);

            end
            %             source = BC{ii};
            q_i{ii,1} = [source ; q_tot_temp((ii-1)*temp+1:ii*temp,:)];

        end



    end

    q_tot = zeros(length(nodes),length(tt));
    for ii = 1:length(RHS1)

        q_tot = q_tot + q_i{ii};

    end

else % Without boundary condition

    for ii = 1:length(RHS1)

        q_i{ii,1} = q_tot((ii-1)*length(nodes)+1:(ii)*length(nodes),:);

    end


    q_tot = zeros(length(nodes),length(tt));
    for ii = 1:length(RHS1)

        q_tot = q_tot + q_i{ii};

    end


end

end

