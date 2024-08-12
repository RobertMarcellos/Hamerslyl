% This case plots results for multivolume discrete GDE (Case 4)
close all
clc
clearvars -except discretization save_format
load([pwd,'/Time_evolutions/MC_discrete_GDE_test_case.mat'])
load([pwd,'/Time_evolutions/MC_discrete_GDE_test_case_evolution.mat'],'PHI','N','v_middle','Vmin','tt')

save_loc = [pwd,'/figs/'];

% In code run invidually, save format have to be uncommented and chosen
% save_format = 'png';

h = fig('width',22,'height',26,'fontsize',13,'Font','Helvetica');

% Choosing discretization
index = find(discretizations == 50);
time_index(1) = find(round(tt,3) == 1);
time_index(2) = find(round(tt,3) == 2);

% Collecting data for plotting
plot_data = [{MC_discrete_FEM_solution(index,1).q_tot},{MC_discrete_FEM_solution(index,1).q1},{MC_discrete_FEM_solution(index,1).q2};
    {MC_discrete_PGFEM_solution(index,1).q_tot},{MC_discrete_PGFEM_solution(index,1).q1},{MC_discrete_PGFEM_solution(index,1).q2};
    {MC_discrete_sec_solution(index,1).q_tot},{MC_discrete_sec_solution(index,1).q1},{MC_discrete_sec_solution(index,1).q2}];

xvector = [{MC_discrete_FEM_solution(index,1).volume},{MC_discrete_PGFEM_solution(index,1).volume},...
    {MC_discrete_sec_solution(index,1).volume}];

plot_colors = [{'blue'},{'red'},{1/256*[102, 51, 0]}];

anal_data = [{PHI(1).phi+PHI(2).phi},{PHI(1).phi},{PHI(2).phi}];

methods = [{'FEM'},{'PGFEM'},{'Sectional'}];
volumes = [{'q_{tot}'},{'q_1'},{'q_2'}];
%% Plotting subplot for temporal evolutions
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

% Plotting errors
subplot(4,3,10)
plot(tt(1:end-1),MC_discrete_FEM_solution(index,1).q_tot_error,'b','linewidth',1.5)
yticks([0,5,10,15,20])
hold on
plot(tt(1:end-1),MC_discrete_PGFEM_solution(index,1).q_tot_error,'r','linewidth',1.5)
plot(tt(1:end-1),MC_discrete_sec_solution(index,1).q_tot_error,'Color',plot_colors{3},'linewidth',1.5)
grid on
xticks([0,0.5,1,1.5,2])
ylabel('Relative error [%]')

subplot(4,3,11)
plot(tt(1:end-1),MC_discrete_FEM_solution(index).q1_error,'b','linewidth',1.5)
yticks([0,5,10,15,20])
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

% saveas(gcf,[save_loc,'fig8'],save_format)


%% Average relative error figure
for kk = 1:size(MC_discrete_sec_solution,2)
    for ii = 1:length(discretizations)




        avg_error_FEM(ii,1) = MC_discrete_FEM_solution(ii,kk).q_tot_avg_error;
        avg_error_PGFEM(ii,1) = MC_discrete_PGFEM_solution(ii,kk).q_tot_avg_error;
        avg_error_sec(ii,1) = MC_discrete_sec_solution(ii,kk).q_tot_avg_error;

        comp_time_FEM(ii,1) = MC_discrete_FEM_solution(ii,kk).evol_time;
        comp_time_PGFEM(ii,1) = MC_discrete_PGFEM_solution(ii,kk).evol_time;
        comp_time_sec(ii,1) = MC_discrete_sec_solution(ii,kk).evol_time;

    end

    FEM{kk,1} = avg_error_FEM;
    PGFEM{kk,1} = avg_error_PGFEM;
    sec{kk,1} = avg_error_sec;


    FEM{kk,2} = comp_time_FEM;
    PGFEM{kk,2} = comp_time_PGFEM;
    sec{kk,2} = comp_time_sec;
    


    clearvars avg_error_FEM avg_error_PGFEM avg_error_sec comp_time_FEM comp_time_PGFEM comp_time_sec

end

h = fig('width',22,'height',10,'fontsize',12,'Font','Helvetica');
subplot(121)
loglog(discretizations,FEM{1,1},'bo','LineWidth',1.2,'MarkerSize',7)
xlim([0,2000])
xticks([1e2 5e2 1e3])
hold on
loglog(discretizations,PGFEM{1,1},'ro','LineWidth',1.2,'MarkerSize',7)
loglog(discretizations,sec{1,1},'o','Color',plot_colors{3},'LineWidth',1.2,'MarkerSize',7)
loglog(discretizations,FEM{2,1},'bx','LineWidth',1.2)
loglog(discretizations,PGFEM{2,1},'rx','LineWidth',1.2)
loglog(discretizations,sec{2,1},'x','Color',plot_colors{3},'LineWidth',1.2)
xlabel([{'Number of bins/elements'};{'(a)'}])
ylabel('Average relative error [%]')
lg = legend('FEM, 0.01h','PGFEM, 0.01h','Sec. 0.01h','FEM, 0.04h','PGFEM, 0.04h','Sec. 0.04h','location','northeast');
lg.FontSize = 10;
lg.NumColumns = 2;
legend boxoff
lg.Position = [0.130201783013033 0.751388870630833 0.336538455162484 0.175496684004929];

subplot(122)
loglog(FEM{1,2},FEM{1,1},'bo--','LineWidth',1.2,'MarkerSize',7)
xticks([1e-1 1e0 1e1 1e2]);
hold on
loglog(PGFEM{1,2},PGFEM{1,1},'ro--','LineWidth',1.2,'MarkerSize',7)
loglog(sec{1,2},sec{1,1},'o--','Color',plot_colors{3},'LineWidth',1.2,'MarkerSize',7)
loglog(FEM{2,2},FEM{2,1},'bx--','LineWidth',1.2)
% hold on
loglog(PGFEM{2,2},PGFEM{2,1},'rx--','LineWidth',1.2)
loglog(sec{2,2},sec{2,1},'x--','Color',plot_colors{3},'LineWidth',1.2)
xlabel([{'Computation time [s]'};{'(b)'}])
ylabel('Average relative error [%]')

% saveas(gcf,[save_loc,'fig9'],save_format)