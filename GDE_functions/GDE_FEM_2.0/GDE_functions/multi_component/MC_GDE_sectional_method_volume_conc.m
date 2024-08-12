function [phi_tot,q_tot,PHI,q,v_middle,v_width,evol_time] = MC_GDE_sectional_method_volume_conc(num_of_bins,v_interval,GR_MC,coag_kernel,removal,source,n_initial,type,t_max,time_step,BC)
%MC_GDE_sectional_method_volume_conc calculates the volume concentrations
%for different substances. In addition, it calculates the volume
%distributions for the time evolution
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2022
%
%
% INPUTS:
% num_of_bins       :   Number of bin in the v_interval
% v_interval        :   Minimun and maximun values for the size interval of
%                       volume of particles as [a,b] where 10^a and 10^b
% GR_MC             :   Growth rate for the different substances as a cell
%                       array containing function for the growth rate for
%                       each cells (If empty, condensation is neglected)
% coag_kernel       :   oagulation kernel function desribing the
%                       coagulation of the particles (If empty, coagulation
%                       is neglected)
% removal           :   Function for the removal term
% source            :   Particle sources into the bins as a cell array
%                       containing vector for source terms
% n_initial         :   Initial distributions for volume concentration or
%                       for volume distribution as a cell array containing
%                       functions
% type              :   Type of initial distribution in order to calculate
%                       volume concentrations for each bin
% t_max             :   End time for the time evolution
% time_step         :   Time step for the time evolution
% BC                :   Possible boundary conditions
%
% OUTPUTS:
% phi_tot           :   total volume concentration of the particles in each
%                       bin
% q_tot             :   total volume distribution
% PHI               :   Individual volume concentrations
% q                 :   Individual volume distributions
% v_middle          :   Middle points for each bin
% v_width           :   Width of the bins
% evol_time         :   Time to calculate temporal evolution with the C-N


% Matrices for the multicomponent sectional time evolution
[v_middle,v_width,C1,C2,Rem,X,Beta] = multicomponent_GDE_sectional_matrices_volume(v_interval,num_of_bins,GR_MC,removal,coag_kernel);

if ~isempty(removal)

    C2 = C2-Rem;

end

% Forming volume concentration from different initial distributions
if strcmp(type,'number')
    for ii = 1:length(n_initial)
        %         N(:,1) = N(:,1) + n_initial(v_middle);
        int_dist = n_initial{ii};

        if isa(int_dist,'function_handle')
            q(ii).q = v_middle.*int_dist(v_middle);
            PHI(ii).phi =  v_width.*q(ii).q;
        else
            q(ii).q = v_middle.*int_dist;
            PHI(ii).phi =  v_width.*int_dist;
        end
    
    end
elseif strcmp(type,'distribution')

    for ii = 1:length(n_initial)
        %         N(:,1) = N(:,1) + v_width*n_initial(v_middle);
        int_dist = n_initial{ii};
        if isa(int_dist,'function_handle')
            q(ii).q = v_width.*int_dist(v_middle);
            PHI(ii).phi =  v_width.*q(ii).q;
        else
            q(ii).q = v_width.*int_dist;
            PHI(ii).phi =  v_width.*int_dist;
        end
        
    
    end

elseif strcmp(type,'volume')
    for ii = 1:length(n_initial)
        %         N(:,1) = N(:,1) + v_width*n_initial(v_middle);
        int_dist = n_initial{ii};
        if isa(int_dist,'function_handle')
            q(ii).q = int_dist(v_middle);
            PHI(ii).phi =  v_width.*int_dist(v_middle);
        else
            q(ii).q = int_dist;
            PHI(ii).phi =  v_width.*int_dist;
        end
    end
end

% Initializing time evolution

time_vec = 0:time_step:t_max;
time_vec = time_vec(:);

% Initializing matrices for the time evolution
if ~isempty(C2)
    G1 = zeros(size(C2));
    % G2 = zeros(size(C2)*length(n_initial));
    for ii = 1:length(n_initial)

        for jj = 1:length(n_initial)
            G1((ii-1)*size(C2,1)+1:(ii*size(C2,1)),(jj-1)*size(C2,1)+1:(jj*size(C2,1))) = C1{ii};

        end
        G2((ii-1)*size(C2,1)+1:(ii*size(C2,1)),(ii-1)*size(C2,1)+1:(ii*size(C2,1))) = C2;

        PHI_tot((ii-1)*length(v_middle)+1:ii*length(v_middle),1) = PHI(ii).phi;

    end



    if nargin > 10

        for ii = length(GR_MC)-1:-1:0

            %             A(length(v_middle)*ii+1,:) = [];
            G1(length(v_middle)*ii+1,:) = [];
            G2(length(v_middle)*ii+1,:) = [];

        end

        K_b = zeros(size(G1,1),1);
        for ii = length(GR_MC)-1:-1:0

            K_b = K_b + G1(:,length(v_middle)*ii+1) + G2(:,length(v_middle)*ii+1);
            G1(:,length(v_middle)*ii+1) = [];
            G2(:,length(v_middle)*ii+1) = [];


        end

        if ~isempty(coag_kernel)
            %             X{1} = [];
            Beta(1,:) = [];
            Beta(:,1) = [];
        end

    end
    % If condensation is only phenomena affecting the time evolution
    S1 = [eye(size(G2))-(time_step/2)*(G1+G2)];
    S2 = [eye(size(G2))+(time_step/2)*(G1+G2)];
else

    for ii = 1:length(n_initial)

        PHI_tot((ii-1)*length(v_middle)+1:ii*length(v_middle),1) = PHI(ii).phi;

    end
end



% Forward integrating through the time evolution

% NOTE: usually boundary condition for sectional method is NOT used
if nargin > 10
    Beta(1,:) = [];
    Beta(:,1) = [];
    if ~exist('K_b','var')

        K_b = zeros(length(n_initial)*(length(v_middle)-1),1);
    
    end
    tic
    for ii = 1:length(time_vec)-1

        % Forming coagulation matrices
        if ~isempty(coag_kernel)


            F = [];
            R = [];

            % Formation matrix F creation
            PHI_sum = zeros(length(v_middle)-1,1);
            for kk = 1:length(n_initial)

                PHI_sum = PHI_sum + PHI_tot((kk-1)*length(v_middle)+2:kk*length(v_middle),end);

            end
%             F = zeros(length(v_middle)-1);
%             for h = 1:length(v_middle)-1
% 
% 
% 
%                 temp = X{h+1};
%                 temp(1,:) = [];
%                 temp(:,1) = [];
%                 F(h,:) = (PHI_sum./v_middle(2:end))'*(Beta.*temp);
% 
% 
%             end


%             R = (PHI_sum).*Beta./v_middle(2:end)';

            for h = 1:length(v_middle)-1

                temp = X{h+1};
                temp(1,:) = [];
                temp(:,1) = [];
                F(h,:) = ((Beta.*temp)*(PHI_sum./v_middle(2:end)))';


            end

            for kk = 1:length(n_initial)

                R{kk} = (PHI_tot((kk-1)*length(v_middle(2:end))+1:kk*length(v_middle(2:end)),end)./v_middle(2:end)').*Beta;

            end


            COAG = zeros(size(F)*length(n_initial));
            for h = 1:length(n_initial)

                COAG((h-1)*length(v_middle(2:end))+1:(h*length(v_middle(2:end))),(h-1)*length(v_middle(2:end))+1:(h*length(v_middle(2:end)))) = F-R{h};
                for kk = 1:length(n_initial)

                    if kk ~= h
                        COAG((h-1)*length(v_middle(2:end))+1:(h*length(v_middle(2:end))),(kk-1)*length(v_middle(2:end))+1:(kk*length(v_middle(2:end)))) = -R{h};

                    end
                end


            end
%             COAG = zeros((length((F))-1)*length(n_initial));
%             for h = 1:length(n_initial)
% 
%                 COAG((h-1)*(length(v_middle)-1)+1:(h*(length(v_middle)-1)),(h-1)*(length(v_middle)-1)+1:(h*(length(v_middle)-1))) = F-R;
% 
% 
%             end


            if ~isempty(C1) % Case for both condensation and coagulation

                S1 = [eye(size(G2))-(time_step/2)*(G1+G2+COAG)];
                S2 = [eye(size(G2))+(time_step/2)*(G1+G2+COAG)];

            else % For Pure Coagulation case

                S1 = [eye(size(COAG))-(time_step/2)*(COAG)];
                S2 = [eye(size(COAG))+(time_step/2)*(COAG)];

            end
            %             S_BC2 = (K_b*BC(time_vec(ii+1)));
            %             S_BC1 = (K_b*BC(time_vec(ii)));
            BC_temp_sum = 0;
            for kk = 1:length(BC)
                BC_temp = BC{kk};
                BC_temp_sum = BC_temp_sum+BC_temp(ii+1);
                S_BC2 = (K_b*BC_temp(ii+1));
                S_BC1 = (K_b*BC_temp(ii));
            end


            PHI_temp = PHI_tot(:,end);
            for jj = length(n_initial)-1:-1:0


                PHI_temp((length(v_middle)-1)*jj+1,:) = [];


            end
            PHI_temp(:,2) = (S1\S2*(PHI_temp))+S1\(S_BC1+S_BC2);


            temp = size(PHI_temp,1)/length(GR_MC);
            for jj = 1:length(GR_MC)



                % Zero condition
                PHI_temp(PHI_temp < 0) = 0;

                PHI_tot((jj-1)*length(v_middle)+1,ii+1) = BC_temp_sum;
                PHI_tot((jj-1)*length(v_middle)+2:(jj)*length(v_middle),ii+1) = PHI_temp((jj-1)*temp+1:jj*temp,end);

            end


            % Placing individuals evolutions for the substances into to the
            % structs
            for h = 1:length(n_initial)

                PHI(h).phi(:,ii+1) = PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);
                q(h).q(:,ii+1) = v_width.^(-1).*PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);


            end


        else

            % Without coagulation S1 and S2 have to be only calculated once

            S_BC2 = (K_b*BC(time_vec(ii+1)));
            S_BC1 = (K_b*BC(time_vec(ii)));


            PHI_temp = PHI_tot(:,end);
            for jj = length(GR_MC)-1:-1:0


                PHI_temp((length(v_middle)-1)*jj+1,:) = [];


            end
            PHI_temp(:,ii+1) = (S1\S2*(PHI_temp))+S1\(S_BC1+S_BC2);


            temp = size(PHI_temp,1)/length(GR_MC);
            for jj = 1:length(GR_MC)


                % Zero condition
                PHI_temp(PHI_temp < 0) = 0;

                PHI_tot((jj-1)*length(v_middle)+1,ii+1) = BC(time_vec(ii+1));
                PHI_tot((jj-1)*length(v_middle)+2:(jj)*length(v_middle),ii+1) = PHI_temp((jj-1)*temp+1:jj*temp,end);

            end

            % Placing individuals evolutions for the substances into to the
            % structs
            for h = 1:length(n_initial)

                PHI(h).phi(:,ii+1) = PHI_tot((h-1)*size(C2,1)+1:(h*size(C2,1)),end);
                q(h).q(:,ii+1) = v_width.^(-1).*PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);


            end


        end


    end
    evol_time = toc;

else

    % Forward integrating through the time evolution
    tic
    for ii = 1:length(time_vec)-1

        % Forming coagulation matrices
        if ~isempty(coag_kernel)

            F = [];
            R = [];
            % Formation matrix F creation
            PHI_sum = zeros(size(v_middle));
            for kk = 1:length(n_initial)

                PHI_sum = PHI_sum + PHI_tot((kk-1)*length(v_middle)+1:kk*length(v_middle),end);

            end
            F = zeros(length(v_middle));
            for h = 1:length(v_middle)


                F(h,:) = ((Beta.*X{h})*(PHI_sum./v_middle))';


            end

            for kk = 1:length(n_initial)

                R{kk} = (PHI_tot((kk-1)*length(v_middle)+1:kk*length(v_middle),end)./v_middle').*Beta;

            end

            COAG = zeros(size(F)*length(n_initial));
            for h = 1:length(n_initial)

                COAG((h-1)*length(v_middle)+1:(h*length(v_middle)),(h-1)*length(v_middle)+1:(h*length(v_middle))) = F-R{h};
                for kk = 1:length(n_initial)

                    if kk ~= h
                        COAG((h-1)*length(v_middle)+1:(h*length(v_middle)),(kk-1)*length(v_middle)+1:(kk*length(v_middle))) = -R{h};

                    end
                end


            end


            if ~isempty(C1) % Case for both condensation and coagulation

                S1 = [eye(size(G2))-(time_step/2)*(G1+G2+COAG)];
                S2 = [eye(size(G2))+(time_step/2)*(G1+G2+COAG)];

            else % For Pure Coagulation case

                S1 = [eye(size(COAG))-(time_step/2)*(COAG)];
                S2 = [eye(size(COAG))+(time_step/2)*(COAG)];

            end

            if ~isempty(source)
                So = zeros(length(source)*length(v_middle),1);
                for jj = 1:length(source)

                    s = source{jj};
                    So((jj-1)*length(v_middle)+1,1) = v_width(1)*s(time_vec(ii))./v_middle(1);
                end

                PHI_tot(:,ii+1) = (S1\S2*(PHI_tot(:,ii)))+(S1\(time_step*So));

            else



                PHI_tot(:,ii+1) = (S1\S2*(PHI_tot(:,ii)));
            end
            % Zero condition
            PHI_tot(PHI_tot < 0) = 0;

            % Placing individuals evolutions for the substances into to the
            % structs
            for h = 1:length(n_initial)

                PHI(h).phi(:,ii+1) = PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);
                q(h).q(:,ii+1) = v_width.^(-1).*PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);


            end


        else

            % Without coagulation S1 and S2 have to be only calculated once

            if ~isempty(source)
                So = zeros(length(source)*length(v_middle),1);
                for jj = 1:length(source)

                    s = source{jj};
                    So((jj-1)*length(v_middle)+1,1) = s(time_vec(ii));
                end

                PHI_tot(:,ii+1) = (S1\S2*(PHI_tot(:,ii)))+time_step*(S1\(So));
            else



                PHI_tot(:,ii+1) = (S1\S2*(PHI_tot(:,ii)));
            end
            % Zero condition
            PHI_tot(PHI_tot < 0) = 0;

            % Placing individuals evolutions for the substances into to the
            % structs
            for h = 1:length(n_initial)

                PHI(h).phi(:,ii+1) = PHI_tot((h-1)*size(C2,1)+1:(h*size(C2,1)),end);
                q(h).q(:,ii+1) = v_width.^(-1).*PHI_tot((h-1)*length(v_middle)+1:(h*length(v_middle)),end);


            end


        end


    end
    evol_time = toc;

end

PHI_sum = zeros(size(q(1).q));
for jj = 1:length(n_initial)


    PHI_sum = PHI_sum + PHI(jj).phi;



end

phi_tot = PHI_sum;
q_tot = v_width.^(-1).*phi_tot;

