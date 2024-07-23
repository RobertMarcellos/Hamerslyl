% This code is used to generate FEM, PGFEM and sectional method solutions
% for the multicomponent GDE test case. These methods are
% compared to the discrete multivolume GDE
clearvars -except discretizations time_steps create_disc_coag create_disc_GDE time_steps_disc
close all
f = waitbar(0,'Hold your horses...');
addpath([pwd,'\Time_evolutions'])

%% CALCULATION TIME EVOLUTIONS

for kk = 1:2
    for jj = 1:length(discretizations)


        clearvars -except discretizations time_steps f jj time_steps_disc MC_discrete_FEM_solution MC_discrete_PGFEM_solution MC_discrete_sec_solution kk create_disc_coag create_disc_GDE
        load('MC_discrete_GDE_test_case_evolution.mat')

        delta_t = time_steps_disc(kk);

%         if kk == 2
% 
%             delta_t = 2;
% 
%         end
        
        % Petrov-Galerkin parameter
        eps = 0.2;
        

        % Coagulation kernel
        fun = @(v,w) 3600*10^(-muunnin)*(vakio*(1./w+1./v).^(1/2).*(v.^(1/3)+w.^(1/3)).^2);
        
        % Initial distribution
        n_initial = @(v)  N1*exp(-0.5*((v-mu1)./sigma1).^2)+N2*exp(-0.5*((v-mu2)./sigma2).^2);

        % Single component GDE
        g = logspace(log10(2*10^(Vmin)),Vmax,discretizations(jj))';

        %% Multicomponent
        initial_dist = [{aa.*g.^(1).*n_initial(g)};{bb.*g.^(1).*n_initial(g)}];

        % Individual growth rates
        GR_MC = [{@(v) 3600*10^(-muunnin)*a*N0.*10^Vmin*vakio*v.^(-1).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2};...
            {@(v) 3600*10^(-muunnin)*bbb*N0.*10^Vmin*vakio*v.^(-1).*(1./(10^Vmin)+1./v).^(1/2).*(v.^(1/3)+(10^Vmin).^(1/3)).^2}];

        [A,G1,G2,B,C,~,A_PG,G1_PG,G2_PG,B_PG,C_PG] = ...
            multicomponent_GDE_FE_matrix_creator(g,GR_MC,[],fun,eps);

        % Boundary conditions
        B1 = 10^-Vmin*PHI(1).phi(1,:)';
        B2 = 10^-Vmin*PHI(2).phi(1,:)';
        dB1 = gradient(B1,delta_t);
        dB2 = gradient(B2,delta_t);
        BC = [{@(ii) B1(ii)},{@(ii) B2(ii)}];
        dBC = [{@(ii) dB1(ii)},{@(ii) dB2(ii)}];
     
        tic
        [q_FEM_tot,q_FEM,nodes,tt] = ...
            multicomponent_GDE_CrankNicolson(g,initial_dist,delta_t,tmax,...
            A,G1,G2,B,C,[],BC,dBC);
        MC_discrete_FEM_solution(jj,kk).evol_time = toc;
        MC_discrete_FEM_solution(jj,kk).q_tot = q_FEM_tot;
        MC_discrete_FEM_solution(jj,kk).volume = g;
        MC_discrete_FEM_solution(jj,kk).q1 = q_FEM{1};
        MC_discrete_FEM_solution(jj,kk).q2 = q_FEM{2};
        MC_discrete_FEM_solution(jj,kk).time_vec = tt;
        MC_discrete_FEM_solution(jj,kk).delta_t = delta_t;



        tic
        [q_PGFEM_tot,q_PGFEM] = ...
            multicomponent_GDE_CrankNicolson(g,initial_dist,delta_t,...
            tmax,A_PG,G1_PG,G2_PG,B_PG,C_PG,[],BC,dBC);
        MC_discrete_PGFEM_solution(jj,kk).evol_time = toc;
        MC_discrete_PGFEM_solution(jj,kk).q_tot = q_PGFEM_tot;
        MC_discrete_PGFEM_solution(jj,kk).volume = g;
        MC_discrete_PGFEM_solution(jj,kk).q1 = q_PGFEM{1};
        MC_discrete_PGFEM_solution(jj,kk).q2 = q_PGFEM{2};
        MC_discrete_PGFEM_solution(jj,kk).time_vec = tt;
        MC_discrete_PGFEM_solution(jj,kk).delta_t = delta_t;
        MC_discrete_PGFEM_solution(jj,kk).epsilon = eps;


        % Sectional method
        initial_dist_sec = [{@(g) aa.*g.^(1).*n_initial(g)};{@(g) bb.*g.^(1).*n_initial(g)}];

        fun_sec = @(v,w) fun(v,w);
        tic
        [phi_tot_sec,q_tot_sec,PHI_sec,q_sec,v_sec,v_width] = MC_GDE_sectional_method_volume_conc(...
            discretizations(jj),[log10(min(d_disc)),Vmax],GR_MC,fun,[],[],initial_dist_sec,'volume',tmax,delta_t);
        MC_discrete_sec_solution(jj,kk).evol_time = toc;
        MC_discrete_sec_solution(jj,kk).q_tot = q_tot_sec;
        MC_discrete_sec_solution(jj,kk).volume = v_sec;
        MC_discrete_sec_solution(jj,kk).volume_width = v_width;
        MC_discrete_sec_solution(jj,kk).q1 = q_sec(1).q;
        MC_discrete_sec_solution(jj,kk).q2 = q_sec(2).q;
        MC_discrete_sec_solution(jj,kk).time_vec = tt;
        MC_discrete_sec_solution(jj,kk).delta_t = delta_t;

        if kk == 1

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


            d_edges = logspace(log10(min(d_disc)),Vmax,discretizations(jj))';
            [MC_discrete_sec_solution(jj,kk).q1_error,MC_discrete_sec_solution(jj,kk).q1_avg_error] = ...
                Error_estimator(log10(v_sec),q_sec(1).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,1:end-1),d_edges,1);
            [MC_discrete_sec_solution(jj,kk).q2_error,MC_discrete_sec_solution(jj,kk).q2_avg_error] = ...
                Error_estimator(log10(v_sec),q_sec(2).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,1:end-1),d_edges,1);
            [MC_discrete_sec_solution(jj,kk).q_tot_error,MC_discrete_sec_solution(jj,kk).q_tot_avg_error] = ...
                Error_estimator(log10(v_sec),q_tot_sec,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,1:end-1)+PHI(2).phi(:,1:end-1)),d_edges,1);

        elseif kk == 2

            tt_disc = round(tt_disc(:),3);
            index = find(mod(tt_disc,delta_t) == 0);

            [MC_discrete_FEM_solution(jj,kk).q1_error,MC_discrete_FEM_solution(jj,kk).q1_avg_error] = ...
                Error_estimator(log10(g),q_FEM{1},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,index));
            [MC_discrete_FEM_solution(jj,kk).q2_error,MC_discrete_FEM_solution(jj,kk).q2_avg_error] = ...
                Error_estimator(log10(g),q_FEM{2},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,index));
            [MC_discrete_FEM_solution(jj,kk).q_tot_error,MC_discrete_FEM_solution(jj,kk).q_tot_avg_error] = ...
                Error_estimator(log10(g),q_FEM_tot,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,index)+PHI(2).phi(:,index)));

            [MC_discrete_PGFEM_solution(jj,kk).q1_error,MC_discrete_PGFEM_solution(jj,kk).q1_avg_error] = ...
                Error_estimator(log10(g),q_PGFEM{1},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,index));
            [MC_discrete_PGFEM_solution(jj,kk).q2_error,MC_discrete_PGFEM_solution(jj,kk).q2_avg_error] = ...
                Error_estimator(log10(g),q_PGFEM{2},log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,index));
            [MC_discrete_PGFEM_solution(jj,kk).q_tot_error,MC_discrete_PGFEM_solution(jj,kk).q_tot_avg_error] = ...
                Error_estimator(log10(g),q_PGFEM_tot,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,index)+PHI(2).phi(:,index)));


            d_edges = logspace(log10(min(d_disc)),Vmax,discretizations(jj))';
            [MC_discrete_sec_solution(jj,kk).q1_error,MC_discrete_sec_solution(jj,kk).q1_avg_error] = ...
                Error_estimator(log10(v_sec),q_sec(1).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(1).phi(:,index),d_edges,1);
            [MC_discrete_sec_solution(jj,kk).q2_error,MC_discrete_sec_solution(jj,kk).q2_avg_error] = ...
                Error_estimator(log10(v_sec),q_sec(2).q,log10(v_middle),diff(v_edges(2,:))^(-1)*PHI(2).phi(:,index),d_edges,1);
            [MC_discrete_sec_solution(jj,kk).q_tot_error,MC_discrete_sec_solution(jj,kk).q_tot_avg_error] = ...
                Error_estimator(log10(v_sec),q_tot_sec,log10(v_middle),diff(v_edges(2,:))^(-1)*(PHI(1).phi(:,index)+PHI(2).phi(:,index)),d_edges,1);



        end
        waitbar(((kk-1)*length(discretizations)+jj)/(2*length(discretizations)+1),...
            f,['Time step: ',num2str(kk),'/',num2str(2),', Discretization: ',num2str(jj),'/',num2str(length(discretizations))])




    end
end

waitbar(1,f,'Finishing');


clearvars -except MC_discrete_sec_solution MC_discrete_FEM_solution MC_discrete_PGFEM_solution discretizations time_steps f time_steps_disc
save([pwd,'\Time_evolutions\MC_discrete_GDE_test_case.mat'],...
    'MC_discrete_sec_solution','MC_discrete_FEM_solution','MC_discrete_PGFEM_solution',...
    'discretizations')

close(f)







