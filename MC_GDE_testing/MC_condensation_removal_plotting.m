% Code for plotting results for test case 2
close all
clearvars -except discretizations save_format
clc

% Load premade data
load([pwd,'/Time_evolutions/MC_condensation_removal_test_case.mat'])
save_loc = [pwd,'/figs/'];
% save_format = 'epsc';
load('viridis_cmap.mat')
set(0,'DefaultFigureColormap',viridis_cmap);

%% FEM/PGFEM comparison
% DATA collection
plot_disc = [125]; % Decision for discretization
time_disc = 0.5; % Decision for time step

% If upper decisions not found, an error will occur
index2 = find(time_disc == time_steps);
index = find(plot_disc == discretizations);
time_vec = cond_rem_analytical_solution(1,index2).time_vec;

% Collecting data for the color plots
v_dense = logspace(log10(min(min(cond_rem_analytical_solution(1,index2).volume))),...
    log10(max(max(cond_rem_analytical_solution(1,index2).volume))),4000)';
for ii = 1:length(time_vec)

    q1_anal(:,ii) = interp1(cond_rem_analytical_solution(1,index2).volume(:,ii),...
        cond_rem_analytical_solution(1,index2).q1(:,ii),...
        v_dense,'linear');
    q2_anal(:,ii) = interp1(cond_rem_analytical_solution(1,index2).volume(:,ii),...
        cond_rem_analytical_solution(1,index2).q2(:,ii),...
        v_dense,'linear');
    q_tot_anal(:,ii) = interp1(cond_rem_analytical_solution(1,index2).volume(:,ii),...
        cond_rem_analytical_solution(1,index2).q_tot(:,ii),...
        v_dense,'linear');
    q1_anal(isnan(q1_anal)) = 0;
    q2_anal(isnan(q2_anal)) = 0;
    q_tot_anal(isnan(q_tot_anal)) = 0;

    q1_FEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_FEM_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_FEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_FEM_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q_tot_FEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_FEM_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');

    q1_PGFEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_PGFEM_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_PGFEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_PGFEM_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q_tot_PGFEM_inter(:,ii) = interp1(cond_rem_FEM_solution(index,index2).volume,...
        abs(cond_rem_PGFEM_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');

    q1_sec_inter(:,ii) = interp1(cond_rem_sec_solution(index,index2).volume,...
        abs(cond_rem_sec_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_sec_inter(:,ii) = interp1(cond_rem_sec_solution(index,index2).volume,...
        abs(cond_rem_sec_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q_tot_sec_inter(:,ii) = interp1(cond_rem_sec_solution(index,index2).volume,...
        abs(cond_rem_sec_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');


end
q1_FEM_inter(isnan(q1_FEM_inter)) = 0;
q2_FEM_inter(isnan(q2_FEM_inter)) = 0;
q_tot_FEM_inter(isnan(q_tot_FEM_inter)) = 0;

q1_PGFEM_inter(isnan(q1_PGFEM_inter)) = 0;
q2_PGFEM_inter(isnan(q2_PGFEM_inter)) = 0;
q_tot_PGFEM_inter(isnan(q_tot_PGFEM_inter)) = 0;

q1_sec_inter(isnan(q1_sec_inter)) = 0;
q2_sec_inter(isnan(q2_sec_inter)) = 0;
q_tot_sec_inter(isnan(q_tot_sec_inter)) = 0;

% Collecting data to cell array
plot_data = [{q_tot_anal},{q1_anal},{q2_anal};
    {q_tot_FEM_inter},{q1_FEM_inter},{q2_FEM_inter};...
    {q_tot_PGFEM_inter},{q1_PGFEM_inter},{q2_PGFEM_inter};...
    {q_tot_sec_inter},{q1_sec_inter},{q2_sec_inter}];

% Collecting average relativ errors and computation times for evolutions
for jj = 1:length(discretizations)

    for ii = 1:length(time_steps)

        avg_error_q_tot_FEM(jj,ii) = cond_rem_FEM_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_FEM(jj,ii) = cond_rem_FEM_solution(jj,ii).q1_avg_error;
        avg_error_q2_FEM(jj,ii) = cond_rem_FEM_solution(jj,ii).q2_avg_error;

        avg_error_q_tot_PGFEM(jj,ii) = cond_rem_PGFEM_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_PGFEM(jj,ii) = cond_rem_PGFEM_solution(jj,ii).q1_avg_error;
        avg_error_q2_PGFEM(jj,ii) = cond_rem_PGFEM_solution(jj,ii).q2_avg_error;

        avg_error_q_tot_sec(jj,ii) = cond_rem_sec_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_sec(jj,ii) = cond_rem_sec_solution(jj,ii).q1_avg_error;
        avg_error_q2_sec(jj,ii) = cond_rem_sec_solution(jj,ii).q2_avg_error;

        q_tot_FEM_comp_time(jj,ii) = cond_rem_FEM_solution(jj,ii).evol_time;
        q_tot_PGFEM_comp_time(jj,ii) = cond_rem_PGFEM_solution(jj,ii).evol_time;
        q_tot_sec_comp_time(jj,ii) = cond_rem_sec_solution(jj,ii).evol_time;


    end


end

[d_grid,t_grid] = meshgrid(v_dense,time_vec);
imgrid.d_grid = d_grid;
imgrid.t_grid = t_grid;

%% PLOTTING

h = fig('width',22,'height',26,'fontsize',12,'Font','Helvetica');

for ii = 1:size(plot_data,1)+1
    for jj = 1:size(plot_data,2)
        subplot(5,3,(ii-1)*size(plot_data,2)+jj)



        if ii == 5

            if jj == 1

                plot(time_vec,cond_rem_FEM_solution(index,index2).q_tot_error,'linewidth',1.5);
                xticks([0 max(time_vec)/2 max(time_vec)])
                xlim([0,max(time_vec)])
                hold on
                plot(time_vec,cond_rem_PGFEM_solution(index,index2).q_tot_error,'linewidth',1.5);
                plot(time_vec,cond_rem_sec_solution(index,index2).q_tot_error,'Color',1/256*[102, 51, 0],'linewidth',1.5);
                ylabel('Relative error [%]','FontSize',13)

                ylim([0,50])
                lg = legend('FEM','PGFEM','Sec.');
                legend boxoff
                lg.FontSize = 10;
                nates = get(gca,"Position");
                lg_nates = lg.Position;
                lg.Position = [nates(1);nates(2)+nates(4)-lg_nates(4);lg_nates(3);lg_nates(4)];

            elseif jj == 2
                plot(time_vec,cond_rem_FEM_solution(index,index2).q1_error,'linewidth',1.5);
                xticks([0 max(time_vec)/2 max(time_vec)])
                xlim([0,max(time_vec)])
                hold on
                plot(time_vec,cond_rem_PGFEM_solution(index,index2).q1_error,'linewidth',1.5);
                plot(time_vec,cond_rem_sec_solution(index,index2).q1_error,'Color',1/256*[102, 51, 0],'linewidth',1.5);
                ylim([0,50])
                %                 ylabel('Relative error [%]')
                lg = legend('FEM','PGFEM','Sec.');
                legend boxoff
                lg.FontSize = 10;
                nates = get(gca,"Position");
                lg_nates = lg.Position;
                lg.Position = [nates(1);nates(2)+nates(4)-lg_nates(4);lg_nates(3);lg_nates(4)];

            elseif jj == 3

                plot(time_vec,cond_rem_FEM_solution(index,index2).q2_error,'linewidth',1.5);
                xticks([0 max(time_vec)/2 max(time_vec)])
                xlim([0,max(time_vec)])
                hold on
                plot(time_vec,cond_rem_PGFEM_solution(index,index2).q2_error,'linewidth',1.5);
                plot(time_vec,cond_rem_sec_solution(index,index2).q2_error,'Color',1/256*[102, 51, 0],'linewidth',1.5);
                ylim([0,50])
                %                 ylabel('Relative error [%]')
                lg = legend('FEM','PGFEM','Sec.');
                legend boxoff
                lg.FontSize = 10;
                nates = get(gca,"Position");
                lg_nates = lg.Position;
                lg.Position = [nates(1);nates(2)+nates(4)-lg_nates(4);lg_nates(3);lg_nates(4)];

            elseif jj == 4

                plot(time_vec,cond_rem_FEM_solution(index,index2).q3_error,'linewidth',1.5);
                xticks([0 max(time_vec)/2 max(time_vec)])
                xlim([0,max(time_vec)])
                hold on
                plot(time_vec,cond_rem_PGFEM_solution(index,index2).q3_error,'linewidth',1.5);
                plot(time_vec,cond_rem_sec_solution(index,index2).q3_error,'Color',1/256*[102, 51, 0],'linewidth',1.5);
                ylim([0,100])
                %                 ylabel('Relative error [%]')
                lg = legend('FEM','PGFEM','Sec.');
                legend boxoff
                lg.FontSize = 7;
                nates = get(gca,"Position");
                lg_nates = lg.Position;
                lg.Position = [nates(1);nates(2)+nates(4)-lg_nates(4);lg_nates(3);lg_nates(4)];

            end

        else
            % Labels for the axis
            labels.ylab = '';
            labels.xlab = '';
            labels.yscale = 'Log';
            labels.ytick = [1e-2 1e-1 1e0 1e1];
            labels.xtick = [0 max(time_vec)/2 max(time_vec)];
            labels.clab = '';
            labels.size = 11;
            labels.title = [];
            fig2.figno = 1;
            fig2.position = [];
            fig2.subplot = [];
            fig2.clf = 0;

            % Data
            ImagesForScaling = [plot_data{ii,jj}];
            c = PlotParticleDensityEvolution(plot_data{ii,jj},imgrid,ImagesForScaling,labels,fig2);
            clim([min(min(q_tot_anal)),max(max(q_tot_anal))])
            colorbar off
            ylim([1e-2,1e1])

            if ii == 1

                if jj == 1
                    title('q_{tot}')
                    ylabel('Analytical','FontSize',13)
                elseif jj == 2
                    title('q_1')
                elseif jj == 3
                    title('q_2')
                else
                    title('q_3')
                end

            elseif ii == 2 && jj == 1

                ylabel('FEM','FontSize',13)

            elseif ii == 3 && jj == 1

                ylabel('PGFEM','FontSize',13)

            elseif ii == 4 && jj == 1

                ylabel('Sectional','FontSize',13)



            end
            if jj == 3 && ii == 4

                c = colorbar('Position',...
                    [0.920028409090909 0.279847756410257 0.0141287878787877 0.644024970981045]);
                c.Label.String = 'Volume concentration [d N/d log(v_p)]';
                c.Label.FontSize = 14;
                clim([min(min(q_tot_anal)),max(max(q_tot_anal))])

            end
        end





    end
end
annotation(h,'textbox',...
    [0.0520321969696968 0.516955128205136 0.364606060606061 0.0386698717948725],...
    'String',{'Particle volume [\mum^3]'},...
    'Rotation',90,...
    'LineStyle','none',...
    'FontSize',13,...
    'FitBoxToText','off');

annotation(h,'textbox',...
    [0.414611742424242 0.0427403853432496 0.205052077580723 0.0345993582464947],...
    'String','Evolution time [h]',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',13,...
    'FitBoxToText','off');
saveas(gcf,[save_loc,'fig4'],save_format)
%%
h = fig('width',22,'height',10,'fontsize',11,'Font','Helvetica');
subplot(121)
semilogy(discretizations,avg_error_q_tot_FEM(:,1),'bo','LineWidth',1.2,'MarkerSize',7)
ylim([0,100])
hold on
loglog(discretizations,avg_error_q_tot_PGFEM(:,1),'ro','LineWidth',1.2,'MarkerSize',7)
loglog(discretizations,avg_error_q_tot_sec(:,1),'o','Color',1/256*[102, 51, 0],'LineWidth',1.2,'MarkerSize',7)
semilogy(discretizations,avg_error_q_tot_FEM(:,2),'bx','LineWidth',1.2)
loglog(discretizations,avg_error_q_tot_PGFEM(:,2),'rx','LineWidth',1.2)
loglog(discretizations,avg_error_q_tot_sec(:,2),'x','Color',1/256*[102, 51, 0],'LineWidth',1.2)
xlabel([{'Number of bins/elements'};{'(a)'}])
ylabel('Average relative error [%]')
lg = legend('FEM, 0.5h','PGFEM, 0.5h','Sec. 0.5h','FEM, 2h','PGFEM, 2h','Sec. 2h','location','northeast');
lg.FontSize = 9.5;
lg.NumColumns = 2;
legend boxoff
lg.Position = [0.173489224272314 0.762633662297499 0.295673071764983 0.175496684004929];

subplot(122)
loglog(q_tot_FEM_comp_time(:,1),avg_error_q_tot_FEM(:,1),'bo--','LineWidth',1.2,'MarkerSize',7)
ylim([0,100])
xticks([1e-1 1e0 1e1 1e2 1e3])
hold on
loglog(q_tot_PGFEM_comp_time(:,1),avg_error_q_tot_PGFEM(:,1),'ro--','LineWidth',1.2,'MarkerSize',7)
loglog(q_tot_sec_comp_time(:,1),avg_error_q_tot_sec(:,1),'o--','Color',1/256*[102, 51, 0],'LineWidth',1.2,'MarkerSize',7)

loglog(q_tot_FEM_comp_time(:,2),avg_error_q_tot_FEM(:,2),'bx--','LineWidth',1.2)
loglog(q_tot_PGFEM_comp_time(:,2),avg_error_q_tot_PGFEM(:,2),'rx--','LineWidth',1.2)
loglog(q_tot_sec_comp_time(:,2),avg_error_q_tot_sec(:,2),'x--','Color',1/256*[102, 51, 0],'LineWidth',1.2)
xlabel([{'Computation time [s]'};{'(b)'}])
ylabel('Average relative error [%]')
saveas(gcf,[save_loc,'fig5'],save_format)