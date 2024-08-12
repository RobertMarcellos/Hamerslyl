% This code plots results for case 3 (multivolume discrete GDE)
close all
clc
clearvars -except discretization save_format

% Loading the data
load([pwd,'/Time_evolutions/MC_discrete_coagulation_test_case.mat'])
load([pwd,'/Time_evolutions/MC_discrete_GDE_coagulation_only.mat'],'PHI','N','v_middle','Vmin','tt')

save_loc = [pwd,'/figs/'];

% If code is run independently save format must be included
% save_format = 'epsc';

h = fig('width',22,'height',26,'fontsize',12,'Font','Helvetica');

% Choose discretization
index = find(discretizations == 50);

% If not found, an error will occurr
time_index(1) = find(round(tt,2) == 1.0);
time_index(2) = find(round(tt,2) == 2.0);

plot_data = [{MC_discrete_FEM_solution(index).q_tot},{MC_discrete_FEM_solution(index).q1},{MC_discrete_FEM_solution(index).q2};
    {MC_discrete_PGFEM_solution(index).q_tot},{MC_discrete_PGFEM_solution(index).q1},{MC_discrete_PGFEM_solution(index).q2};
    {MC_discrete_sec_solution(index).q_tot},{MC_discrete_sec_solution(index).q1},{MC_discrete_sec_solution(index).q2}];

xvector = [{MC_discrete_FEM_solution(index).volume},{MC_discrete_PGFEM_solution(index).volume},...
    {MC_discrete_sec_solution(index).volume}];

plot_colors = [{'blue'},{'red'},{1/256*[102, 51, 0]}];

anal_data = [{PHI(1).phi+PHI(2).phi},{PHI(1).phi},{PHI(2).phi}];

methods = [{'FEM'},{'PGFEM'},{'Sectional'}];
volumes = [{'q_{tot}'},{'q_1'},{'q_2'}];
%% Plotting temporal evolutions at 48h and 06h
for jj = 1:length(methods)
    for ii = 1:length(volumes)


        subplot(4,3,(jj-1)*length(methods)+ii)
        semilogx(v_middle,10^(-Vmin)*anal_data{ii}(:,1),'k--','LineWidth',1.3)
        xticks([1e-8,1e-7,1e-6])
        hold on
        for kk = 1:length(time_index)

            semilogx(v_middle,10^(-Vmin)*anal_data{ii}(:,time_index(kk)),'k--','LineWidth',2.5)
            if jj == 3
                semilogx(xvector{jj},plot_data{jj,ii}(:,time_index(kk)),'Color',plot_colors{jj},'linewidth',1.3);
            else
                semilogx(xvector{jj},plot_data{jj,ii}(:,time_index(kk)),plot_colors{jj},'linewidth',1.3);
            end



        end

        if jj == 1

            title(volumes{ii})
        end

        if jj == 3 & ii == 2

            xlabel('Particle volume [\mum^3]')

        end
        if jj == 2 & ii == 1

            mod_label = [{'Volume concentration [dN/d log(v_p)]'};{methods{jj}}];
            ylabel(mod_label)

        elseif ii == 1

            ylabel(methods(jj))

        end

    end
end


%% Plotting errors for temporal evolutions
subplot(4,3,10)
plot(tt(1:end-1),MC_discrete_FEM_solution(index).q_tot_error,'b','linewidth',1.5)
% yticks([0,5,10,15,20])
hold on
plot(tt(1:end-1),MC_discrete_PGFEM_solution(index).q_tot_error,'r','linewidth',1.5)
plot(tt(1:end-1),MC_discrete_sec_solution(index).q_tot_error,'Color',plot_colors{3},'linewidth',1.5)
grid on
xticks([0,0.5,1,1.5,2])
ylabel('Relative error [%]')

subplot(4,3,11)
plot(tt(1:end-1),MC_discrete_FEM_solution(index).q1_error,'b','linewidth',1.5)
% yticks([0,5,10,15,20])
hold on
plot(tt(1:end-1),MC_discrete_PGFEM_solution(index).q1_error,'r','linewidth',1.5)
plot(tt(1:end-1),MC_discrete_sec_solution(index).q1_error,'Color',plot_colors{3},'linewidth',1.5)
grid on

xticks([0,0.5,1,1.5,2])
xlabel('Evolution time [h]')

subplot(4,3,12)
plot(tt(1:end-1),MC_discrete_FEM_solution(index).q2_error,'b','linewidth',1.5)
% yticks([0,5,10,15,20])
hold on
plot(tt(1:end-1),MC_discrete_PGFEM_solution(index).q2_error,'r','linewidth',1.5)
plot(tt(1:end-1),MC_discrete_sec_solution(index).q2_error,'Color',plot_colors{3},'linewidth',1.5)
grid on

xticks([0,0.5,1,1.5,2])
saveas(gcf,[save_loc,'fig6'],save_format)


%% Average relative error figure

for ii = 1:length(discretizations)

    avg_error_FEM(ii,1) = MC_discrete_FEM_solution(ii).q_tot_avg_error;
    avg_error_PGFEM(ii,1) = MC_discrete_PGFEM_solution(ii).q_tot_avg_error;
    avg_error_sec(ii,1) = MC_discrete_sec_solution(ii).q_tot_avg_error;

    avg_error_FEM_q1(ii,1) = MC_discrete_FEM_solution(ii).q1_avg_error;
    avg_error_PGFEM_q1(ii,1) = MC_discrete_PGFEM_solution(ii).q1_avg_error;
    avg_error_sec_q1(ii,1) = MC_discrete_sec_solution(ii).q1_avg_error;


    avg_error_FEM_q2(ii,1) = MC_discrete_FEM_solution(ii).q2_avg_error;
    avg_error_PGFEM_q2(ii,1) = MC_discrete_PGFEM_solution(ii).q2_avg_error;
    avg_error_sec_q2(ii,1) = MC_discrete_sec_solution(ii).q2_avg_error;

    comp_time_FEM(ii,1) = MC_discrete_FEM_solution(ii).evol_time;
    comp_time_PGFEM(ii,1) = MC_discrete_PGFEM_solution(ii).evol_time;
    comp_time_sec(ii,1) = MC_discrete_sec_solution(ii).evol_time;

end

h = fig('width',22,'height',10,'fontsize',12,'Font','Helvetica');
subplot(121)
loglog(discretizations,avg_error_FEM(:,1),'bo','LineWidth',1.2)
hold on
loglog(discretizations,avg_error_PGFEM(:,1),'rx','LineWidth',1.2)
loglog(discretizations,avg_error_sec(:,1),'*','Color',plot_colors{3},'LineWidth',1.2)
xlabel([{'Number of bins/elements'};{'(a)'}])
ylabel('Average relative error [%]')
legend('FEM','PGFEM','Sec.','location','northeast')
legend boxoff
% title('q_{tot}')

% subplot(222)
% loglog(discretizations,avg_error_FEM_q1(:,1),'bo-','LineWidth',1.5)
% hold on
% loglog(discretizations,avg_error_PGFEM_q1(:,1),'rx-','LineWidth',1.5)
% loglog(discretizations,avg_error_sec_q1(:,1),'*-','Color',plot_colors{3},'LineWidth',1.5)
% xlabel([{'Number of bins/elements'};{'(b)'}])
% ylabel('Average relative error [%]')
% legend('FEM','PGFEM','Sec.','location','northeast')
% legend boxoff
% title('q_1')
% 
% subplot(223)
% loglog(discretizations,avg_error_FEM_q2(:,1),'bo-','LineWidth',1.5)
% hold on
% loglog(discretizations,avg_error_PGFEM_q2(:,1),'rx-','LineWidth',1.5)
% loglog(discretizations,avg_error_sec_q2(:,1),'*-','Color',plot_colors{3},'LineWidth',1.5)
% xlabel([{'Number of bins/elements'};{'(c)'}])
% ylabel('Average relative error [%]')
% legend('FEM','PGFEM','Sec.','location','northeast')
% legend boxoff
% title('q_2')

subplot(122)
loglog(comp_time_FEM(:,1),avg_error_FEM(:,1),'bo--','LineWidth',1.2)
hold on
loglog(comp_time_PGFEM(:,1),avg_error_PGFEM(:,1),'rx--','LineWidth',1.2)
loglog(comp_time_sec(:,1),avg_error_sec(:,1),'*--','Color',plot_colors{3},'LineWidth',1.2)
xlabel([{'Computation time [s]'};{'(b)'}])
ylabel('Average relative error [%]')
legend('FEM','PGFEM','Sec.','location','northeast')
legend boxoff

saveas(gcf,[save_loc,'fig7'],save_format)