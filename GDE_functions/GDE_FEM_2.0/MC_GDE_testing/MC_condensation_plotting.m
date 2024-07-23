% Code for plotting results for test case 1 (Only condensation affecting
% aerosol number distribution)
close all
clearvars -except save_format discretizations
clc

% Load premade time evolutions
load([pwd,'/Time_evolutions/MC_condensation_test_case.mat'])
save_loc = [pwd,'/figs/'];

load('viridis_cmap.mat')
set(0,'DefaultFigureColormap',viridis_cmap);

% Uncomment if this code is run independently
% save_format = 'epsc';

%% FEM/PGFEM comparison

% DATA collection
plot_disc = [125]; % Choose discretization for plotting
time_disc = 0.5; % Choose time discretization

% Decisions are look from the data (if they do not exist, this will give on
% error)
index2 = find(time_disc == time_steps);
index = find(plot_disc == discretizations);
time_vec = cond_analytical_solution(1,index2).time_vec;

% Colorplot initializations
v_dense = logspace(log10(min(min(cond_analytical_solution(1,index2).volume))),...
    log10(max(max(cond_analytical_solution(1,index2).volume))),4000)';
for ii = 1:length(time_vec)

    q1_anal(:,ii) = interp1(cond_analytical_solution(1,index2).volume(:,ii),...
        cond_analytical_solution(1,index2).q1(:,ii),...
        v_dense,'linear');
    q2_anal(:,ii) = interp1(cond_analytical_solution(1,index2).volume(:,ii),...
        cond_analytical_solution(1,index2).q2(:,ii),...
        v_dense,'linear');
    q3_anal(:,ii) = interp1(cond_analytical_solution(1,index2).volume(:,ii),...
        cond_analytical_solution(1,index2).q3(:,ii),...
        v_dense,'linear');
    q_tot_anal(:,ii) = interp1(cond_analytical_solution(1,index2).volume(:,ii),...
        cond_analytical_solution(1,index2).q_tot(:,ii),...
        v_dense,'linear');
    q1_anal(isnan(q1_anal)) = 0;
    q2_anal(isnan(q2_anal)) = 0;
    q3_anal(isnan(q3_anal)) = 0;
    q_tot_anal(isnan(q_tot_anal)) = 0;

    q1_FEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_FEM_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_FEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_FEM_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q3_FEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_FEM_solution(index,index2).q3(:,ii)),...
        v_dense,'linear');
    q_tot_FEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_FEM_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');

    q1_PGFEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_PGFEM_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_PGFEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_PGFEM_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q3_PGFEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_PGFEM_solution(index,index2).q3(:,ii)),...
        v_dense,'linear');
    q_tot_PGFEM_inter(:,ii) = ...
        interp1(cond_FEM_solution(index,index2).volume,...
        abs(cond_PGFEM_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');

    q1_sec_inter(:,ii) = ...
        interp1(cond_sec_solution(index,index2).volume,...
        abs(cond_sec_solution(index,index2).q1(:,ii)),...
        v_dense,'linear');
    q2_sec_inter(:,ii) = ...
        interp1(cond_sec_solution(index,index2).volume,...
        abs(cond_sec_solution(index,index2).q2(:,ii)),...
        v_dense,'linear');
    q3_sec_inter(:,ii) = ...
        interp1(cond_sec_solution(index,index2).volume,...
        abs(cond_sec_solution(index,index2).q3(:,ii)),...
        v_dense,'linear');
    q_tot_sec_inter(:,ii) = ...
        interp1(cond_sec_solution(index,index2).volume,...
        abs(cond_sec_solution(index,index2).q_tot(:,ii)),...
        v_dense,'linear');


end
q1_FEM_inter(isnan(q1_FEM_inter)) = 0;
q2_FEM_inter(isnan(q2_FEM_inter)) = 0;
q3_FEM_inter(isnan(q3_FEM_inter)) = 0;
q_tot_FEM_inter(isnan(q_tot_FEM_inter)) = 0;

q1_PGFEM_inter(isnan(q1_PGFEM_inter)) = 0;
q2_PGFEM_inter(isnan(q2_PGFEM_inter)) = 0;
q3_PGFEM_inter(isnan(q3_PGFEM_inter)) = 0;
q_tot_PGFEM_inter(isnan(q_tot_PGFEM_inter)) = 0;

q1_sec_inter(isnan(q1_sec_inter)) = 0;
q2_sec_inter(isnan(q2_sec_inter)) = 0;
q3_sec_inter(isnan(q3_sec_inter)) = 0;
q_tot_sec_inter(isnan(q_tot_sec_inter)) = 0;

% Collecting plottable data to a cell array
plot_data = [{q_tot_anal},{q1_anal},{q2_anal},{q3_anal};
    {q_tot_FEM_inter},{q1_FEM_inter},{q2_FEM_inter},{q3_FEM_inter};...
    {q_tot_PGFEM_inter},{q1_PGFEM_inter},{q2_PGFEM_inter},{q3_PGFEM_inter};...
    {q_tot_sec_inter},{q1_sec_inter},{q2_sec_inter},{q3_sec_inter}];

colorbar_lims = [min(min(q_tot_anal)),max(max(q_tot_anal));
    min(min(q1_anal)),max(max(q1_anal));
    min(min(q2_anal)),max(max(q2_anal));
    min(min(q3_anal)),max(max(q3_anal))];

concentration = [{'q_{tot}'},{'q_1'},{'q_2'},{'q_{3}'}];
methods = [{'Analytical'},{'FEM'},{'PGFEM'},{'Sectional'}];

% Collecting average relative errors and computation times for evolutions
for jj = 1:length(discretizations)

    for ii = 1:length(time_steps)

        avg_error_q_tot_FEM(jj,ii) = cond_FEM_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_FEM(jj,ii) = cond_FEM_solution(jj,ii).q1_avg_error;
        avg_error_q2_FEM(jj,ii) = cond_FEM_solution(jj,ii).q2_avg_error;
        avg_error_q3_FEM(jj,ii) = cond_FEM_solution(jj,ii).q3_avg_error;

        avg_error_q_tot_PGFEM(jj,ii) = cond_PGFEM_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_PGFEM(jj,ii) = cond_PGFEM_solution(jj,ii).q1_avg_error;
        avg_error_q2_PGFEM(jj,ii) = cond_PGFEM_solution(jj,ii).q2_avg_error;
        avg_error_q3_PGFEM(jj,ii) = cond_PGFEM_solution(jj,ii).q3_avg_error;

        avg_error_q_tot_sec(jj,ii) = cond_sec_solution(jj,ii).q_tot_avg_error;
        avg_error_q1_sec(jj,ii) = cond_sec_solution(jj,ii).q1_avg_error;
        avg_error_q2_sec(jj,ii) = cond_sec_solution(jj,ii).q2_avg_error;
        avg_error_q3_sec(jj,ii) = cond_sec_solution(jj,ii).q3_avg_error;

        q_tot_FEM_comp_time(jj,ii) = cond_FEM_solution(jj,ii).evol_time;
        q_tot_PGFEM_comp_time(jj,ii) = cond_PGFEM_solution(jj,ii).evol_time;
        q_tot_sec_comp_time(jj,ii) = cond_sec_solution(jj,ii).evol_time;


    end


end

[d_grid,t_grid] = meshgrid(v_dense,time_vec);
imgrid.d_grid = d_grid;
imgrid.t_grid = t_grid;

%% PLOTTING

h = fig('width',22,'height',26,'fontsize',12,'Font','Helvetica');

for ii = 1:size(plot_data,1)
    for jj = 1:size(plot_data,2)
        h = subplot(4,4,(ii-1)*size(plot_data,2)+jj);

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
        % Plot are done with PlotParticleDensityEvolutions function
        c = PlotParticleDensityEvolution(plot_data{ii,jj},imgrid,ImagesForScaling,labels,fig2);
        ylim([1e-2,1e1]);
        clim([colorbar_lims(jj,:)])
        if ii == 4

            
            colorbar off
            pause(0.1)
            loc_pic = get(h,'Position');
            pause(0.1)
            c = colorbar('SouthOutside');
            h.Position = loc_pic;
            h.Position(2) = loc_pic(2)+0.02; 
            c_pos = c.Position;
            pause(0.1)
            c.Position = [loc_pic(1),loc_pic(2)-0.035,loc_pic(3),0.007];
        else
            colorbar off
            pause(0.1)
            loc_pic = get(h,'Position');
            h.Position(2) = loc_pic(2)+0.02; 
        end
        if ii == 1
            title(concentration{jj});

        end

        if jj == 1

            ylabel(methods{ii},'FontSize',13)
        end
        xlabel('Evolution time [h]')

    end


end
annotation('textbox',...
    [0 0 1 0.0376522427708962],...
    'String',{'Volume concentration [d N/d log(v_p)]'},...
    'HorizontalAlignment','center',...
    'FontSize',13,...
    'EdgeColor','none');

annotation('textbox',...
    [0.0578229174014228 0.398910257411078 0.280819120777235 0.0437580118196993],...
    'String','Particle volume [\mum^3]',...
    'Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',13,...
    'EdgeColor','none');
% annotation('textbox',...
%     [0.989479166666667 0.327676282051286 0.461818181818182 0.0284935897435889],...
%     'String','Volume concentration [d N/d log(v_p)]',...
%     'Rotation',90,...
%     'HorizontalAlignment','center',...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

saveas(gcf,[save_loc,'fig2'],save_format)
%% Plotting error figures
h = fig('width',22,'height',26,'fontsize',12,'Font','Helvetica');
subplot(321)
plot(cond_FEM_solution(index,index2).time_vec,cond_FEM_solution(index,index2).q_tot_error,'b','LineWidth',2)
yticks([0,20,40,60,80,100])
hold on
plot(cond_PGFEM_solution(index,index2).time_vec,cond_PGFEM_solution(index,index2).q_tot_error,'r','LineWidth',2)
plot(cond_sec_solution(index,index2).time_vec,cond_sec_solution(index,index2).q_tot_error,'Color',1/256*[102, 51, 0],'LineWidth',2)
xlabel([{'Evolution time [h]'};{'(a)'}])
ylabel('Relative error [%]')
ylim([0,100])
title('q_{tot}')
lg = legend('FEM','PGFEM','Sectional','Location','northwest');
lg.FontSize = 11;
legend boxoff

subplot(322)
plot(cond_FEM_solution(index,index2).time_vec,cond_FEM_solution(index,index2).q1_error,'b','LineWidth',2)
yticks([0,20,40,60,80,100])
hold on
plot(cond_PGFEM_solution(index,index2).time_vec,cond_PGFEM_solution(index,index2).q1_error,'r','LineWidth',2)
plot(cond_sec_solution(index,index2).time_vec,cond_sec_solution(index,index2).q1_error,'Color',1/256*[102, 51, 0],'LineWidth',2)
xlabel([{'Evolution time [h]'};{'(b)'}])
ylabel('Relative error [%]')
title('q_1')
ylim([0,100])
lg = legend('FEM','PGFEM','Sectional','Location','northwest');
lg.FontSize = 11;
legend boxoff

subplot(323)
plot(cond_FEM_solution(index,index2).time_vec,cond_FEM_solution(index,index2).q2_error,'b','LineWidth',2)
yticks([0,20,40,60,80,100])
hold on
plot(cond_PGFEM_solution(index,index2).time_vec,cond_PGFEM_solution(index,index2).q2_error,'r','LineWidth',2)
plot(cond_sec_solution(index,index2).time_vec,cond_sec_solution(index,index2).q2_error,'Color',1/256*[102, 51, 0],'LineWidth',2)
xlabel([{'Evolution time [h]'};{'(c)'}])
ylabel('Relative error [%]')
title('q_2')
ylim([0,100])
lg = legend('FEM','PGFEM','Sectional','Location','northwest');
lg.FontSize = 11;
legend boxoff

subplot(324)
plot(cond_FEM_solution(index,index2).time_vec,cond_FEM_solution(index,index2).q3_error,'b','LineWidth',2)
yticks([0,20,40,60,80,100])
hold on
plot(cond_PGFEM_solution(index,index2).time_vec,cond_PGFEM_solution(index,index2).q3_error,'r','LineWidth',2)
plot(cond_sec_solution(index,index2).time_vec,cond_sec_solution(index,index2).q3_error,'Color',1/256*[102, 51, 0],'LineWidth',2)
xlabel([{'Evolution time [h]'};{'(d)'}])
ylabel('Relative error [%]')
title('q_3')
ylim([0,100])
lg = legend('FEM','PGFEM','Sectional','Location','northwest');
lg.FontSize = 11;
legend boxoff

subplot(325)
semilogy(discretizations,avg_error_q_tot_FEM(:,index2),'bo','LineWidth',1.2,'MarkerSize',7)
ylim([0,100])
hold on
loglog(discretizations,avg_error_q_tot_PGFEM(:,index2),'ro','LineWidth',1.2,'MarkerSize',7)
loglog(discretizations,avg_error_q_tot_sec(:,index2),'o','Color',1/256*[102, 51, 0],'LineWidth',1.2,'MarkerSize',7)
semilogy(discretizations,avg_error_q_tot_FEM(:,2),'bx','LineWidth',1.2)
loglog(discretizations,avg_error_q_tot_PGFEM(:,2),'rx','LineWidth',1.2)
loglog(discretizations,avg_error_q_tot_sec(:,2),'x','Color',1/256*[102, 51, 0],'LineWidth',1.2)
xlabel([{'Number of bins/elements'};{'(e)'}])
ylabel('Average relative error [%]')
lg = legend('FEM, 0.5h','PGFEM, 0.5h','Sec. 0.5h','FEM, 2h','PGFEM, 2h','Sec. 2h','location','northeast');
lg.FontSize = 9.5;
lg.NumColumns = 2;
legend boxoff
lg.Position = [0.167479608887698 0.266875474236115 0.295673071764983 0.0539165804369161];

subplot(326)
loglog(q_tot_FEM_comp_time(:,index2),avg_error_q_tot_FEM(:,index2),'bo--','LineWidth',1.2,'MarkerSize',7)
ylim([0,100])
xticks([1e-1 1e0 1e1 1e2 1e3])
hold on
loglog(q_tot_PGFEM_comp_time(:,index2),avg_error_q_tot_PGFEM(:,index2),'ro--','LineWidth',1.2,'MarkerSize',7)
loglog(q_tot_sec_comp_time(:,index2),avg_error_q_tot_sec(:,index2),'o--','Color',1/256*[102, 51, 0],'LineWidth',1.2,'MarkerSize',7)

loglog(q_tot_FEM_comp_time(:,2),avg_error_q_tot_FEM(:,2),'bx--','LineWidth',1.2)
loglog(q_tot_PGFEM_comp_time(:,2),avg_error_q_tot_PGFEM(:,2),'rx--','LineWidth',1.2)
loglog(q_tot_sec_comp_time(:,2),avg_error_q_tot_sec(:,2),'x--','Color',1/256*[102, 51, 0],'LineWidth',1.2)
xlabel([{'Computation time [s]'};{'(f)'}])
ylabel('Average relative error [%]')

saveas(gcf,[save_loc,'fig3'],save_format)


%% Example figure of real temporal evolutions
time_vector = cond_analytical_solution(index2).time_vec;
d_disc = cond_analytical_solution(index2).volume;
time_index(1) = find(time_vector == 84);
time_index(2) = find(time_vector == 168);
h = fig('width',22,'height',12,'fontsize',13,'Font','Helvetica');
semilogx(d_disc(:,1),cond_analytical_solution(index2).q_tot(:,1),...
    'k','LineWidth',2.5,'DisplayName','q_{tot}')
xlim([1e-2,1e1]);
hold on
semilogx(d_disc(:,1),cond_analytical_solution(index2).q1(:,1),...
    'r','LineWidth',1.5,'DisplayName','q_1')
semilogx(d_disc(:,1),cond_analytical_solution(index2).q2(:,1),...
    'b','LineWidth',1.5,'DisplayName','q_2')
semilogx(d_disc(:,1),cond_analytical_solution(index2).q3(:,1),...
    'm','LineWidth',1.5,'DisplayName','q_3')

h1 = semilogx(d_disc(:,time_index(1)),...
    cond_analytical_solution(index2).q_tot(:,time_index(1)),...
    'k--','LineWidth',2.5,'DisplayName','q_{tot}');
h2 = semilogx(d_disc(:,time_index(1)),...
    cond_analytical_solution(index2).q1(:,time_index(1)),...
    'r--','LineWidth',1.5,'DisplayName','q_1');
h3 = semilogx(d_disc(:,time_index(1)),...
    cond_analytical_solution(index2).q2(:,time_index(1)),...
    'b--','LineWidth',1.5,'DisplayName','q_2');
h4 = semilogx(d_disc(:,time_index(1)),...
    cond_analytical_solution(index2).q3(:,time_index(1)),...
    'm--','LineWidth',1.5,'DisplayName','q_3');

h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
h4.Annotation.LegendInformation.IconDisplayStyle = 'off';

h1 = semilogx(d_disc(:,time_index(2)),...
    cond_analytical_solution(index2).q_tot(:,time_index(2)),...
    'k:','LineWidth',2.5,'DisplayName','q_{tot}');
h2 = semilogx(d_disc(:,time_index(2)),...
    cond_analytical_solution(index2).q1(:,time_index(2)),...
    'r:','LineWidth',1.5,'DisplayName','q_1');
h3 = semilogx(d_disc(:,time_index(2)),...
    cond_analytical_solution(index2).q2(:,time_index(2)),...
    'b:','LineWidth',1.5,'DisplayName','q_2');
h4 = semilogx(d_disc(:,time_index(2)),...
    cond_analytical_solution(index2).q3(:,time_index(2)),...
    'm:','LineWidth',1.5,'DisplayName','q_3');

h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
h4.Annotation.LegendInformation.IconDisplayStyle = 'off';

lg = legend('location','best');
legend boxoff
xlabel('Particle volume [\mum^3]')
ylabel(['Volume concentration [d N/d log(v_p)]'])

annotation('textbox',...
    [0.456501893939394 0.477447916666666 0.109643939393939 0.0683506944444439],...
    'String',{'t = 0h'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.691018939393939 0.678090277777775 0.109643939393939 0.0683506944444441],...
    'String','t = 84h',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
    [0.843755681818182 0.87432291666666 0.109643939393939 0.0683506944444441],...
    'String','t = 168h',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

saveas(gcf,[save_loc,'fig1'],save_format)


