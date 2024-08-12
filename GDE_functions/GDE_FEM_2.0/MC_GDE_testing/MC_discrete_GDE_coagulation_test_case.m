% This code is used to generate FEM, PGFEM and sectional method solutions
% for the multicomponent GDE test case. These methods are
% compared to the discrete multivolume GDE
clearvars -except discretizations time_steps create_disc_coag create_disc_GDE time_steps_disc
close all
f = waitbar(0,'Hold your horses...');
addpath([pwd,'\Time_evolutions'])

%% CALCULATION TIME EVOLUTIONS


for jj = 1:length(discretizations)

    clearvars -except discretizations time_steps time_steps_disc f jj MC_discrete_FEM_solution MC_discrete_PGFEM_solution MC_discrete_sec_solution create_disc_coag create_disc_GDE
    load('MC_discrete_GDE_coagulation_only.mat')

    % Petrov-Galerkin parameter
    eps = 0.2;
    kk = 1;

    % Brownian coagulation kernel
    fun = @(v,w) 3600*10^(-muunnin)*(vakio*(1./w+1./v).^(1/2).*(v.^(1/3)+w.^(1/3)).^2);

    % Initial distribution
    n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2)+N2*exp(-0.5*((v-mu2)./sigma2).^2);

    % Single component GDE
    g = logspace(log10(2*10^(Vmin)),Vmax,discretizations(jj))';

    %% Multicomponent

    % Parameters for components are downloaded from mat-file
    aa = k*g+b;
    aa(aa < y2) = y2;
    bb = 1-aa;

    % Create initial distribution
    initial_dist = [{aa.*g.^(1).*n_initial(g)};{bb.*g.^(1).*n_initial(g)}];

    % Growth rate for components
    GR_MC = [{@(v) 3600*10^(-muunnin)*a*N0.*10^Vmin*vakio*v.^(-1).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2};...
        {@(v) 3600*10^(-muunnin)*bbb*N0.*10^Vmin*vakio*v.^(-1).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2}];

    % Creating FEM matrices
    [A,G1,G2,B,C,~,A_PG,G1_PG,G2_PG,B_PG,C_PG] = ...
        multicomponent_GDE_FE_matrix_creator(g,GR_MC,[],fun,eps);

    B1 = 10^-Vmin*PHI(1).phi(1,:)';
    B2 = 10^-Vmin*PHI(2).phi(1,:)';
    dB1 = gradient(B1,delta_t);
    dB2 = gradient(B2,delta_t);
    BC = [{@(ii) B1(ii)},{@(ii) B2(ii)}];
    dBC = [{@(ii) dB1(ii)},{@(ii) dB2(ii)}];

    tic
    [q_FEM_tot,q_FEM,nodes,tt] = ...
        multicomponent_GDE_CrankNicolson(g,initial_dist,delta_t,tmax,A,...
        {zeros(size(G2));zeros(size(G2))},zeros(size(G2)),B,C,[],BC,dBC);
    MC_discrete_FEM_solution(jj,kk).evol_time = toc;
    MC_discrete_FEM_solution(jj,kk).q_tot = q_FEM_tot;
    MC_discrete_FEM_solution(jj,kk).volume = g;
    MC_discrete_FEM_solution(jj,kk).q1 = q_FEM{1};
    MC_discrete_FEM_solution(jj,kk).q2 = q_FEM{2};
    MC_discrete_FEM_solution(jj,kk).time_vec = tt;
    MC_discrete_FEM_solution(jj,kk).delta_t = delta_t;

    tic
    [q_PGFEM_tot,q_PGFEM] = ...
        multicomponent_GDE_CrankNicolson(g,initial_dist,delta_t,tmax,...
        A_PG,{zeros(size(G2));zeros(size(G2))},zeros(size(G2)),B_PG,...
        C_PG,[],BC,dBC);
    MC_discrete_PGFEM_solution(jj,kk).evol_time = toc;
    MC_discrete_PGFEM_solution(jj,kk).q_tot = q_PGFEM_tot;
    MC_discrete_PGFEM_solution(jj,kk).volume = g;
    MC_discrete_PGFEM_solution(jj,kk).q1 = q_PGFEM{1};
    MC_discrete_PGFEM_solution(jj,kk).q2 = q_PGFEM{2};
    MC_discrete_PGFEM_solution(jj,kk).time_vec = tt;
    MC_discrete_PGFEM_solution(jj,kk).delta_t = delta_t;
    MC_discrete_PGFEM_solution(jj,kk).epsilon = eps;

    % Sectional method
    fit = @(v) k*v+b;
    initial_dist_sec = [{@(g) max(fit(g),y2*ones(size(g))).*g.^(1).*n_initial(g)};{@(g) min(1-fit(g),(1-y2)*ones(size(g))).*g.^(1).*n_initial(g)}];

    fun_sec = @(v,w) fun(v,w);
    tic
    [phi_tot_sec,q_tot_sec,PHI_sec,q_sec,v_sec,v_width] = ...
        MC_GDE_sectional_method_volume_conc(...
        discretizations(jj),[log10(min(v_middle)),Vmax],[],fun,[],[],...
        initial_dist_sec,'volume',tmax,delta_t);
    MC_discrete_sec_solution(jj,kk).evol_time = toc;
    MC_discrete_sec_solution(jj,kk).q_tot = q_tot_sec;
    MC_discrete_sec_solution(jj,kk).volume = v_sec;
    MC_discrete_sec_solution(jj,kk).volume_width = v_width;
    MC_discrete_sec_solution(jj,kk).q1 = q_sec(1).q;
    MC_discrete_sec_solution(jj,kk).q2 = q_sec(2).q;
    MC_discrete_sec_solution(jj,kk).time_vec = tt;
    MC_discrete_sec_solution(jj,kk).delta_t = delta_t;

    [MC_discrete_FEM_solution(jj,kk).q1_error,MC_discrete_FEM_solution(jj,kk).q1_avg_error] = ...
        Error_estimator(log10(g),q_FEM{1},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,1:end-1));
    [MC_discrete_FEM_solution(jj,kk).q2_error,MC_discrete_FEM_solution(jj,kk).q2_avg_error] = ...
        Error_estimator(log10(g),q_FEM{2},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,1:end-1));
    [MC_discrete_FEM_solution(jj,kk).q_tot_error,MC_discrete_FEM_solution(jj,kk).q_tot_avg_error] = ...
        Error_estimator(log10(g),q_FEM_tot,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,1:end-1)+PHI(2).phi(:,1:end-1)));

    [MC_discrete_PGFEM_solution(jj,kk).q1_error,MC_discrete_PGFEM_solution(jj,kk).q1_avg_error] = ...
        Error_estimator(log10(g),q_PGFEM{1},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,1:end-1));
    [MC_discrete_PGFEM_solution(jj,kk).q2_error,MC_discrete_PGFEM_solution(jj,kk).q2_avg_error] = ...
        Error_estimator(log10(g),q_PGFEM{2},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,1:end-1));
    [MC_discrete_PGFEM_solution(jj,kk).q_tot_error,MC_discrete_PGFEM_solution(jj,kk).q_tot_avg_error] = ...
        Error_estimator(log10(g),q_PGFEM_tot,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,1:end-1)+PHI(2).phi(:,1:end-1)));


    d_edges = logspace(log10(min(v_middle)),Vmax,discretizations(jj))';
    [MC_discrete_sec_solution(jj,kk).q1_error,MC_discrete_sec_solution(jj,kk).q1_avg_error] = ...
        Error_estimator(log10(v_sec),q_sec(1).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,1:end-1),d_edges,1);
    [MC_discrete_sec_solution(jj,kk).q2_error,MC_discrete_sec_solution(jj,kk).q2_avg_error] = ...
        Error_estimator(log10(v_sec),q_sec(2).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,1:end-1),d_edges,1);
    [MC_discrete_sec_solution(jj,kk).q_tot_error,MC_discrete_sec_solution(jj,kk).q_tot_avg_error] = ...
        Error_estimator(log10(v_sec),q_tot_sec,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,1:end-1)+PHI(2).phi(:,1:end-1)),d_edges,1);

    waitbar((jj)/(length(discretizations)+1),f,['Discretization: ',num2str(jj),'/',num2str(length(discretizations))])



end

waitbar(1,f,'Finishing');


clearvars -except MC_discrete_sec_solution MC_discrete_FEM_solution MC_discrete_PGFEM_solution discretizations time_steps time_steps_disc f create_disc_GDE
save([pwd,'/Time_evolutions/MC_discrete_coagulation_test_case.mat'],...
    'MC_discrete_sec_solution','MC_discrete_FEM_solution','MC_discrete_PGFEM_solution',...
    'discretizations')

close(f)







