% This code is used to generate FEM, PGFEM and sectional method solutions
% for the multicomponent condensation equation case. These methods are
% compared to the analytical solution presented in
clearvars -except discretizations time_steps create_disc_coag create_disc_GDE time_steps_disc
close all
f = waitbar(0,'Hold your horses...');
%% Inializations for the time evolution calculations

% Size interval (as \mum^3)
Vmin = -3;
Vmax = 3;

% Discretization for the analytical solution
vv = logspace(Vmin,Vmax,3000)';

% Upwinding factor for the PGFEM
epsilon = 0.2;

% Density of particles as g/cm3:
% global rho_p
rho_p = 1;

% Functions used in the variable transformations
diam_to_vol = @(d) pi*d.^(3)/6; % mum --> mum^3 (Spherical particles)
vol_to_diam = @(v) (6*v/pi).^(1/3); % mum^3 --> mum
vol_to_mass = @(v,rho) rho*10^(-6)*v; % mum^3 --> mug
mass_to_vol = @(m,rho) (rho*10^(-6))^(-1)*m; % mug --> mum^3

% Initial volume concentrations
sigma1 = 3e-2;
mu1 = 1e-1;
N1 = 80000;
q_1_initial = @(v) N1*(sigma1*sqrt(2*pi))^(-1)*exp(-0.5*((v-mu1)/sigma1).^2);

sigma2 = 6.8e-2;
mu2 = 2.3e-1;
N2 = 300000;
q_2_initial = @(v) N2*(sigma2*sqrt(2*pi))^(-1)*exp(-0.5*((v-mu2)/sigma2).^2);

q_tot = q_1_initial(vv)+q_2_initial(vv);


% GROWTH RATES for various species
I1 = vol_to_mass(diam_to_vol(3.5),rho_p);
I2 = vol_to_mass(diam_to_vol(4.5),rho_p);

I_tot = I1 + I2;

I = [I_tot;I1;I2];

% Time evolution parameters
t_max = 168;

%% CALCULATION TIME EVOLUTIONS

for kk = 1:length(time_steps)

    clearvars tt A G1 G2 A_PG G1_PG G2_PG m_FEM
    tt = [0:time_steps(kk):t_max]';

    for jj = 1:length(discretizations)



        %% Analytical solution is calculated only once
        if jj == 1

            clearvars Vppt
            % Analytical solution is obtained as function of paricle mass
            % (Analytical solution presented in Katoshevski and Seinfeld
            % 1997)
            m = vol_to_mass(vv,rho_p);
            delta = 2/3;
            m0 = m;
            Vppt(:,1) = vv;
            q0 = q_1_initial(vv)+q_2_initial(vv);
            q01 = q_1_initial(vv);
            q02 = q_2_initial(vv);
            tic
            for ii = 1:length(tt)-1

                % Calculating new mass axis during the evolution
                m_t = [];
                m_t = (m0.^(1-delta)+(1-delta)*I_tot*tt(ii+1)).^(1/(1-delta));
                m(:,ii+1) = m_t;
                Vppt(:,ii+1) = mass_to_vol(m(:,end),rho_p);

            end


            % Applying constant removal rate
            removal = @(m) 3e-3*ones(size(m));
            R = cumtrapz(tt,removal(m)');
            R_tilde = exp(-R');

            for ii = 1:length(tt)-1

                q01(:,ii+1) = R_tilde(:,ii+1).*((m0./m(:,ii+1)).^delta...
                    .*((I1/I_tot)*(m(:,ii+1)./m0-1)+(q01(:,1)./q0(:,1))))...
                    .*q0(:,1);
                q02(:,ii+1) = R_tilde(:,ii+1).*((m0./m(:,ii+1)).^delta...
                    .*((I2/I_tot)*(m(:,ii+1)./m0-1)+(q02(:,1)./q0(:,1))))...
                    .*q0(:,1);

            end

            % Collecting analytical solutions for volume concentrations
            cond_rem_analytical_solution(1,kk).evol_time = toc;
            cond_rem_analytical_solution(1,kk).volume = Vppt;
            cond_rem_analytical_solution(1,kk).q1 = q01;
            cond_rem_analytical_solution(1,kk).q2 = q02;
            cond_rem_analytical_solution(1,kk).q_tot = q01+q02;
            cond_rem_analytical_solution(1,kk).mass = m;
            cond_rem_analytical_solution(1,kk).time_vec = tt;
            cond_rem_analytical_solution(1,kk).delta_t = time_steps(kk);

        end
        waitbar(((kk-1)*length(discretizations)+jj)/(length(time_steps)...
            *length(discretizations)+1),f,...
            ['Time step: ',num2str(kk),'/',num2str(length(time_steps)),...
            ', Discretization: ',num2str(jj),'/',num2str(length(discretizations))])


        %% Estimation method approksimations

        % FEM INITIALIZATIONS
        m_FEM = vol_to_mass(logspace(Vmin,Vmax,discretizations(jj))',rho_p);
        growthrate = [{@(m) I1.*m.^(-1/3)},{@(m) I2.*m.^(-1/3)}];
        [A,G1,G2,~,~,R,A_PG,G1_PG,G2_PG,R_PG] = ...
            multicomponent_GDE_FE_matrix_creator(m_FEM,growthrate,[],[],...
            epsilon,@(m) removal(m));
        init_dist = [{q_1_initial(mass_to_vol(m_FEM,rho_p))};...
            {q_2_initial(mass_to_vol(m_FEM,rho_p))}];

        % Time evolution
        tic
        % removal(m_FEM,rho_p,mass_to_vol,vol_to_diam)

        % NOTE: Removal rate is easier to apply this way if it same for
        % each compound
        [q_tot_FEM,PHI,nodes] = ...
            multicomponent_GDE_CrankNicolson(m_FEM,init_dist,...
            time_steps(kk),t_max,A,G1,G2-A*diag(removal(m_FEM)),...
            [],[],[],[{@(t) 0};{@(t) 0}],[{@(t) 0};{@(t) 0}]);
        time_FEM = toc;
        tic
        [q_tot_PGFEM,PHI_PG] = ...
            multicomponent_GDE_CrankNicolson(m_FEM,init_dist,time_steps(kk),...
            t_max,A_PG,G1_PG,G2_PG-A_PG*diag(removal(m_FEM)),...
            [],[],[],[{@(t) 0};{@(t) 0}],[{@(t) 0};{@(t) 0}]);
        time_PGFEM = toc;
        v_FEM = mass_to_vol(m_FEM,rho_p);

        % Sectional method
        d_edges = vol_to_mass(logspace(Vmin,Vmax,discretizations(jj))',rho_p);
        d_widths = diff(d_edges);
        m_diff = d_edges(1:end-1) + .5*d_widths; % midpoints of the bins
        init_dist = [{q_1_initial(mass_to_vol(m_diff,rho_p))};...
            {q_2_initial(mass_to_vol(m_diff,rho_p))}];

        [phi_tot_sec,q_tot_sec,PHI_sec,q_sec,m_middle,~,time_sec] = ...
            MC_GDE_sectional_method_volume_conc(discretizations(jj)...
            ,[log10(min(m_diff)),log10(max(m_diff))],growthrate,[],...
            @(m) removal(m),[],init_dist,'volume',t_max,time_steps(kk),@(t) 0);
        v_middle = mass_to_vol(m_middle,rho_p);

        % Collecting variables for approximation methods
        cond_rem_FEM_solution(jj,kk).q_tot = q_tot_FEM;
        cond_rem_FEM_solution(jj,kk).volume = v_FEM;
        cond_rem_FEM_solution(jj,kk).mass = m_FEM;
        cond_rem_FEM_solution(jj,kk).q1 = PHI{1};
        cond_rem_FEM_solution(jj,kk).q2 = PHI{2};
        cond_rem_FEM_solution(jj,kk).time_vec = tt;
        cond_rem_FEM_solution(jj,kk).delta_t = time_steps(kk);
        cond_rem_FEM_solution(jj,kk).evol_time = time_FEM;


        cond_rem_PGFEM_solution(jj,kk).q_tot = q_tot_PGFEM;
        cond_rem_PGFEM_solution(jj,kk).volume = v_FEM;
        cond_rem_PGFEM_solution(jj,kk).mass = m_FEM;
        cond_rem_PGFEM_solution(jj,kk).q1 = PHI_PG{1};
        cond_rem_PGFEM_solution(jj,kk).q2 = PHI_PG{2};
        cond_rem_PGFEM_solution(jj,kk).time_vec = tt;
        cond_rem_PGFEM_solution(jj,kk).delta_t = time_steps(kk);
        cond_rem_PGFEM_solution(jj,kk).evol_time = time_PGFEM;

        cond_rem_sec_solution(jj,kk).q_tot = q_tot_sec;
        cond_rem_sec_solution(jj,kk).volume = v_middle;
        cond_rem_sec_solution(jj,kk).mass = m_middle;
        cond_rem_sec_solution(jj,kk).q1 = q_sec(1).q;
        cond_rem_sec_solution(jj,kk).q2 = q_sec(2).q;
        cond_rem_sec_solution(jj,kk).time_vec = tt;
        cond_rem_sec_solution(jj,kk).delta_t = time_steps(kk);
        cond_rem_sec_solution(jj,kk).evol_time = time_sec;


        Error_FEM_q1 = [];
        Error_FEM_q2 = [];
        Error_FEM_q_tot = [];


        Error_PGFEM_q1 = [];
        Error_PGFEM_q2 = [];
        Error_PGFEM_q_tot = [];

        Error_sec_q1 = [];
        Error_sec_q2 = [];
        Error_sec_q_tot = [];


        % Calculating relative errors
        for ii = 1:length(tt)
            v_true = cond_rem_analytical_solution(1,kk).volume(:,ii);

            index_FE1 = find(v_FEM >= min(v_true));
            index_FE2 = find(v_true<= v_FEM(end));
            index1 = find(v_true >= v_FEM(index_FE1(1)) & v_true <= v_true(index_FE2(end)));
            

            index_sec1 = find(v_middle >= min(v_true));
            index_sec2 = find(v_true<= v_middle(end));
            index2 = find(v_true >= v_middle(index_sec1(1)) & v_true <= v_true(index_sec2(end)));

            Vp_temp_FE = v_true(index1);
            q_tot_temp_FE = cond_rem_analytical_solution(1,kk).q_tot(index1,ii);
            q1_temp_FE = cond_rem_analytical_solution(1,kk).q1(index1,ii);
            q2_temp_FE = cond_rem_analytical_solution(1,kk).q2(index1,ii);

            Vp_temp_sec = v_true(index2);
            q_tot_temp_sec = cond_rem_analytical_solution(1,kk).q_tot(index2,ii);
            q1_temp_sec = cond_rem_analytical_solution(1,kk).q1(index2,ii);
            q2_temp_sec = cond_rem_analytical_solution(1,kk).q2(index2,ii);

            q_tot_temp_FE(isnan(q_tot_temp_FE)) = 0;
            q1_temp_FE(isnan(q1_temp_FE)) = 0;
            q2_temp_FE(isnan(q2_temp_FE)) = 0;


            q_tot_temp_sec(isnan(q_tot_temp_sec)) = 0;
            q1_temp_sec(isnan(q1_temp_sec)) = 0;
            q2_temp_sec(isnan(q2_temp_sec)) = 0;
            
            temp1 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_FEM_solution(jj,kk).q1(index_FE1(1):end,ii)...
                ,log10(Vp_temp_FE),q1_temp_FE);
            temp2 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_FEM_solution(jj,kk).q2(index_FE1(1):end,ii),...
                log10(Vp_temp_FE),q2_temp_FE);
            temp7 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_PGFEM_solution(jj,kk).q1(index_FE1(1):end,ii),...
                log10(Vp_temp_FE),q1_temp_FE);
            temp8 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_PGFEM_solution(jj,kk).q2(index_FE1(1):end,ii)...
                ,log10(Vp_temp_FE),q2_temp_FE);

            temp4 = Error_estimator(log10(v_middle(index_sec1(1):end)),...
                cond_rem_sec_solution(jj,kk).q1(index_sec1(1):end,ii),...
                log10(Vp_temp_sec),q1_temp_sec);
            temp5 = Error_estimator(log10(v_middle(index_sec1(1):end)),...
                cond_rem_sec_solution(jj,kk).q2(index_sec1(1):end,ii),...
                log10(Vp_temp_sec),q2_temp_sec);

            temp10 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_FEM_solution(jj,kk).q_tot(index_FE1(1):end,ii),...
                log10(Vp_temp_FE),q_tot_temp_FE);
            temp11 = Error_estimator(log10(v_FEM(index_FE1(1):end)),...
                cond_rem_PGFEM_solution(jj,kk).q_tot(index_FE1(1):end,ii),...
                log10(Vp_temp_FE),q_tot_temp_FE);
            temp12 = Error_estimator(log10(v_middle(index_sec1(1):end)),...
                cond_rem_sec_solution(jj,kk).q_tot(index_sec1(1):end,ii),...
                log10(Vp_temp_sec),q_tot_temp_sec);

            % Collecting errors
            Error_FEM_q1(ii) = temp1;
            Error_FEM_q2(ii) = temp2;
            Error_PGFEM_q1(ii) = temp7;
            Error_PGFEM_q2(ii) = temp8;
            Error_sec_q1(ii) = temp4;
            Error_sec_q2(ii) = temp5;
            Error_FEM_q_tot(ii) = temp10;
            Error_PGFEM_q_tot(ii) = temp11;
            Error_sec_q_tot(ii) = temp12;

        end
        cond_rem_FEM_solution(jj,kk).q1_error =  Error_FEM_q1;
        cond_rem_FEM_solution(jj,kk).q2_error =  Error_FEM_q2;
        cond_rem_FEM_solution(jj,kk).q_tot_error =  Error_FEM_q_tot;
        cond_rem_FEM_solution(jj,kk).q1_avg_error =  mean(Error_FEM_q1);
        cond_rem_FEM_solution(jj,kk).q2_avg_error =  mean(Error_FEM_q2);
        cond_rem_FEM_solution(jj,kk).q_tot_avg_error =  mean(Error_FEM_q_tot);


        cond_rem_PGFEM_solution(jj,kk).q1_error =  Error_PGFEM_q1;
        cond_rem_PGFEM_solution(jj,kk).q2_error =  Error_PGFEM_q2;
        cond_rem_PGFEM_solution(jj,kk).q_tot_error =  Error_PGFEM_q_tot;

        cond_rem_PGFEM_solution(jj,kk).q1_avg_error =  mean(Error_PGFEM_q1);
        cond_rem_PGFEM_solution(jj,kk).q2_avg_error =  mean(Error_PGFEM_q2);
        cond_rem_PGFEM_solution(jj,kk).q_tot_avg_error =  mean(Error_PGFEM_q_tot);

        cond_rem_sec_solution(jj,kk).q1_error =  Error_sec_q1;
        cond_rem_sec_solution(jj,kk).q2_error =  Error_sec_q2;
        cond_rem_sec_solution(jj,kk).q_tot_error =  Error_sec_q_tot;

        cond_rem_sec_solution(jj,kk).q1_avg_error =  mean(Error_sec_q1);
        cond_rem_sec_solution(jj,kk).q2_avg_error =  mean(Error_sec_q2);
        cond_rem_sec_solution(jj,kk).q_tot_avg_error =  mean(Error_sec_q_tot);



    end

end
waitbar(1,f,'Finishing');


clearvars -except cond_rem_analytical_solution cond_rem_sec_solution cond_rem_FEM_solution cond_rem_PGFEM_solution discretizations time_steps f time_steps_disc
save([pwd,'\Time_evolutions\MC_condensation_removal_test_case.mat'],...
    'cond_rem_analytical_solution','cond_rem_sec_solution','cond_rem_FEM_solution',...
    'cond_rem_PGFEM_solution','discretizations','time_steps')

close(f)







