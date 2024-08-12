function [N,PHI,v_middle,v_edges,tt] ...
    = discrete_multivolume_GDE(v_monomer,v_int,initial,type,v_distribution,...
    monomer_num,coag_kernel,time_step,t_max,phenom,removal,source)
%DISCRETE_MULTIVOLUME_GDE calculates the discrete time evolution of the volume
%concentration of different substances based on the GDE. In addition, the
%code can be used to calculate the standard discrete GDE under the
%assumption that particles consist only a single substance.
%
% Teemu Salminen
% University of Eastern Finland
% Department of Technical Physics
% 2022
%
%
%   The codes allows any number of subtances to be included into the model.
%   However, the computational cost increases linearly when more substances
%   are added.
%
%   INPUTS:
%   v_monomer   :   Volume of a monomer (in current version, the volume have
%                   to be same for each substance)
%   v_int       :   Minimun and maximum size of particles as a vector
%                   [Vmin,Vmax]. Note! Vmin have to be a manifold of the
%                   v_monomer
%   initial     :   Initial size distribution n or number concentration N
%   type        :   If initial is size distribution 'distribution'. If number
%                   concentration 'number'
%   v_distribution: How many percent of volume in each size class is of a
%                   substance i. Length can be a single number or length of
%                   v_middle.
%   monomer_num :   Number of monomers for each substance as a vector.
%                   Affects the condensation rate.
%   coag_kernel :   Coagulation kernel function for the condensation and
%                   coagulation.
%   time_step   :   time step for the time evolution calculation
%   t_max       :   The time wen the evolution calculations are ended.
%   phenom      :   set as 'GDE','cond' or 'coag' . Corresponding phenomena
%                   are calculated
%   mono        :   Removal process as a function
%   source      :   Source for new particles
%
%   OUTPUTS
%   N           :   Number concentration calculated from volume
%                   concentrations of different substances
%   PHI         :   Volume concentrations for different substances
%   v_middle    :   Spatial volume vector for the particle volume (center
%                   points of the bins)
%   v_edges     :   Edges for the bins
%   tt          :   Time vector for the time evolution

%% Initializing the size axis and initial ditribution

% Smallest particle size is 2
if v_int(1) == v_monomer
    v_middle = [v_int(1)+v_int(1):v_monomer:v_int(2)]';
else
    v_middle = [v_int(1):v_monomer:v_int(2)]';
end
v_width = v_middle(3)-v_middle(2);
v_edges = [v_middle-0.5*v_width,v_middle+0.5*v_width];

% Creating initial distributions
if strcmp(type,'number')
    N(:,1) = initial(v_middle);
elseif strcmp(type,'distribution')
    N(:,1) = v_width*initial(v_middle);
end

for ii = 1:length(v_distribution)

    if length(v_distribution{ii}) == 1

        PHI(ii).phi = v_distribution{ii}*v_middle.*N;

    else

        temp = v_distribution{ii};
        if ii == 1
            temp(isnan(temp)) = 0;
        elseif ii == 2
            temp(isnan(temp)) = 1;
        end

        qi_temp = temp.*v_middle.*(N./v_width);
        PHI(ii).phi = v_width*qi_temp;

    end

end

%% Condensation and coagulation matrices for multicomponent GDE

if strcmp(phenom,'cond') || strcmp(phenom,'GDE')

    % Condensation matrix (Have to be only formed once)
    C_tot = zeros(length(monomer_num)*length(v_middle),length(monomer_num)*length(v_middle));

    % Applying removal process
    if nargin > 10
        if ~isempty(removal)

            Rem = zeros(length(monomer_num)*length(v_middle),length(monomer_num)*length(v_middle));


        end

    end

    % Creating condensation matrices
    for ii = 1:length(monomer_num)
        for jj = 1:length(monomer_num)

            C_temp = [];
            R_temp = [];

            if ii == jj

                C_temp = diag(-sum(monomer_num).*coag_kernel(v_monomer,v_middle)) + ...
                    diag((sum(monomer_num)+v_monomer*v_middle(1:end-1).^(-1)*monomer_num(ii)).*coag_kernel(v_monomer,v_middle(1:end-1)),-1);
            else

                C_temp = zeros(length(v_middle)) + ...
                    diag((v_monomer*v_middle(1:end-1).^(-1)*monomer_num(ii)).*coag_kernel(v_monomer,v_middle(1:end-1)),-1);

            end

            if nargin > 10
                if ~isempty(removal)

                    if ii == jj
                        R_temp = diag(removal(v_middle));
                    end



                end

            end

            C_tot(length(v_middle)*(ii-1)+1:ii*length(v_middle),...
                length(v_middle)*(jj-1)+1:jj*length(v_middle)) = ...
                C_temp;
            if nargin > 10
                if ~isempty(removal)

                    if ii == jj
                        Rem(length(v_middle)*(ii-1)+1:ii*length(v_middle),...
                            length(v_middle)*(jj-1)+1:jj*length(v_middle)) = ...
                            R_temp;
                    end

                end
            end



        end
    end

    if nargin > 10
        if ~isempty(removal)

            C_tot = C_tot-Rem;

        end
    end
    S = eye(length(C_tot))+(time_step/2)*(C_tot);
    S2 = eye(length(C_tot))-(time_step/2)*(C_tot);


elseif nargin > 10

    if ~isempty(removal)

        for ii = 1:length(monomer_num)
            for jj = 1:length(monomer_num)
                if ii == jj
                    R_temp = diag(removal(v_middle));

                    Rem(length(v_middle)*(ii-1)+1:ii*length(v_middle),...
                        length(v_middle)*(jj-1)+1:jj*length(v_middle)) = ...
                        R_temp;
                end

            end
        end

    end

end

% Creating coagulation matrices
% NOTE: kommentoituna toinen mahdollinen tapa muodostaa
% koagulaatiomatriisit. Lopputuloksena hiukan erilainen matriisiyhtälö,
% joka antaa käytännössä saman tuloksen
if strcmp(phenom,'coag') || strcmp(phenom,'GDE')


    % Coagulation formation and removal matrices (These parts have to be formed only once)
    [~,index1]=min(abs(v_middle-2*v_middle(1)));
    index2 = find(v_middle == v_middle(end)-v_middle(1));

    F_vol = zeros(length(v_middle));
    R_vol = zeros(length(v_middle));

    for ii = 1:index2

        F_vol(index1+ii-1:end,ii) = ...
            coag_kernel(v_middle(ii),v_middle(1:index2+1-ii))./v_middle(ii);
    end

    for ii = 1:length(v_middle)


        R_vol(:,ii) = coag_kernel(v_middle(ii),v_middle)./v_middle(ii);

    end
end

% Initializing the vector for the volume concentrations
phi_temp = [];
% phi_sum = zeros(size(v_middle));
for jj = 1:length(monomer_num)


    phi_temp = [phi_temp ; PHI(jj).phi(:,end)];

end

%% Time evolution calculations for the multicomponent case
tt = [0:time_step:t_max]';
for ii = 1:length(tt)-1

    if strcmp(phenom,'coag') || strcmp(phenom,'GDE')
        % Coagulation matrices depend on the previous values in time
        coag_vol_form = cell(length(monomer_num),1);
        coag_vol_rem = cell(length(monomer_num),1);

        % Values depends from the previous values of all substances
        for kk = 1:length(monomer_num)
            temp = [];
            temp = phi_temp(length(v_middle)*(kk-1)+1:kk*length(v_middle),end);
            temp2 = zeros(size(F_vol));
            for jj = 1:index2

                temp2(index1+jj-1:end,jj) = temp(1:index2+1-jj,1);
            end

            coag_vol_form{kk,1} = temp2.*F_vol;
            coag_vol_rem{kk,1} = -temp.*R_vol;

        end

        % Coagulation matrix is formed from the previously calculated values
        coag_tot = zeros(length(monomer_num)*length(v_middle));
        for kk = 1:length(monomer_num)
            for jj = 1:length(monomer_num)

                coag_tot(length(v_middle)*(kk-1)+1:kk*length(v_middle),...
                    length(v_middle)*(jj-1)+1:jj*length(v_middle)) = ...
                    coag_vol_form{kk,1}+coag_vol_rem{kk,1};

            end
        end

    end

    if strcmp(phenom,'GDE')
        % Crank-Nicolson matrices for case with coagulation and condensation
        S = eye(length(C_tot))+(time_step/2)*(C_tot+coag_tot);
        S2 = eye(length(C_tot))-(time_step/2)*(C_tot+coag_tot);

    elseif strcmp(phenom,'coag')
        % Crank-Nicolson matrices for case with coagulation and condensation

        if nargin > 10
            if ~isempty(removal)

                S = eye(length(coag_tot))+(time_step/2)*(coag_tot-Rem);
                S2 = eye(length(coag_tot))-(time_step/2)*(coag_tot-Rem);

            else

                S = eye(length(coag_tot))+(time_step/2)*(coag_tot);
                S2 = eye(length(coag_tot))-(time_step/2)*(coag_tot);

            end

        else
            S = eye(length(coag_tot))+(time_step/2)*(coag_tot);
            S2 = eye(length(coag_tot))-(time_step/2)*(coag_tot);

        end
    end


    % Calculating volume concentration values for the next time step
    if nargin > 11
        So = zeros(length(source)*length(v_middle),1);
        for jj = 1:length(source)

            s = source{jj};
            So((jj-1)*length(v_middle)+1,1) = s(tt(ii));
        end

        % Note: If source is time dependent this has to be done differently
        phi_temp(:,ii+1) = (S2\S)*phi_temp(:,ii)+S2\(time_step*So);
    else
        phi_temp(:,ii+1) = (S2\S)*phi_temp(:,ii);
        % Euler method for testing
        %         phi_temp(:,ii+1) = phi_temp(:,ii) + ...
        %           time_step*(C_tot+coag_tot)*phi_temp(:,ii);
    end

    % Possible positivity constraint

    %     if sum(find(phi_temp < 0)) > 0
    %         disp('Negative value found')
    %         pause
    %     end
    %         phi_temp(phi_temp < 0) = 0;

    % Calculating the number concentration from volume concentrations
    phi_sum = zeros(size(v_middle));
    for jj = 1:length(monomer_num)


        phi_sum = [phi_sum + phi_temp(length(v_middle)*(jj-1)+1:jj*length(v_middle),end)];

    end
    N(:,ii+1) = phi_sum./v_middle;


end

% Separating the volume concentrations, molecules of different substance
% in some size classes and the current particle compositions
for ii = 1:length(monomer_num)

    temp = phi_temp(length(v_middle)*(ii-1)+1:ii*length(v_middle),:);
    PHI(ii).phi = temp;
    PHI(ii).particle = temp/v_monomer;
    PHI(ii).composition = temp./(v_monomer.*N);


end


end


