% This code plot figures for the chosen time evolutions
close all
clc
% Decision for the test cases (If 1 plot figs for the test case).

% 1. condensation
% 2. coagulation
% 3. GDE analytical
% 4. discrete GDE
% plot_case = [1 1 1 1];
% % Save location for the subplots
% sl_sub = [pwd,'\figs\subplots\'];
% sl_fig = [pwd,'\figs\fig_files\'];

% Choose save format for the subplots
% save_format = 'png';
% save_format = 'epsc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Condensation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_case(1) == 1
    close all
    clearvars -except plot_case sl_sub save_format discretization
    load('condensation_time_evolutions.mat')
    
    if ~exist([pwd,'\figs\figs_condensation'],'dir')
        mkdir([pwd,'\figs\figs_condensation'])
    end
    save_loc = [pwd,'\figs\figs_condensation\'];
    
    % Number of elements/bins for the approximation method
    num_of_disc = discretization(1);
    
    for ii = 2:size(cond_time_evolution,1)
        
        % FE approximation
        cond_hila(ii-1,1) = length(cond_time_evolution{ii,14});
        cond_error_t1(ii-1,1) = cond_time_evolution{ii,7};
        cond_error_t2(ii-1,1) = cond_time_evolution{ii,18};
        
        cond_FEM_t1(ii-1,1) = cond_time_evolution{ii,3};
        
        % Sectional method
        cond_error_t1_diff(ii-1,1) = cond_time_evolution{ii,11};
        cond_error_t2_diff(ii-1,1) = cond_time_evolution{ii,22};
        
        % Implicit Euler approximations
        cond_FEM_ie(ii-1,1) = cond_time_evolution{ii,30};
        cond_PGFEM_ie(ii-1,1) = cond_time_evolution{ii,34};
        
    end
    
    for jj = 1:2 % Standard figures and subplot
        
        if jj == 1
            h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            H = fig('width',20,'height',20','font','Helvetica','fontsize',11)
            subplot(224)
            
        end
        semilogy(cond_hila,cond_error_t1,'bo--','linewidth',1.5)
        hold on
        semilogy(cond_hila,cond_error_t2,'cx--','linewidth',1.5)
        semilogy(cond_hila,cond_error_t1_diff,'rs--','linewidth',1.5)
        semilogy(cond_hila,cond_error_t2_diff,'md--','linewidth',1.5)
        ylabel('Average relative error [%]')
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        lgd = legend('PGFEM, \Deltat = 1h','PGFEM, \Deltat = 0.1h','Sec., \Deltat = 1h',...
            'Sec., \Deltat = 0.1h','location','best');
        lgd.FontSize = 9,5;
        legend boxoff
        set(lgd,'Position',...
            [0.70542401643925 0.358309881040843 0.201058197076674...
            0.091931214408269]);
        
        if jj == 1
            xlabel('Number of bins/elements')
            saveas(gcf,[save_loc,'condensation_average_error.png'])
            h2 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Number of bins/elements'};{'(d)'}]);
            subplot(223)
            
        end
        loglog(cond_hila,cond_error_t1,'bo--','linewidth',1.5)
        hold on
        loglog(cond_hila,cond_FEM_t1,'cx--','linewidth',1.5)
        loglog(cond_hila,cond_PGFEM_ie,'md--','linewidth',1.5)
        loglog(cond_hila,cond_FEM_ie,'rs--','linewidth',1.5)
        ylabel('Average relative error [%]')
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        lgd = legend('C-N PGFEM','C-N FEM','I.E PGFEM','I.E FEM','location','best');
        lgd.FontSize = 9,5;
        legend boxoff
        %         saveas(gcf,[save_loc,'condensation_average_error.png'])
        
        if jj == 1
            xlabel('Number of bins/elements')
            saveas(gcf,[save_loc,'condensation_C-N_IE_comparison.png'])
            h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Number of bins/elements'};{'(c)'}]);
            subplot(221)
            
        end
        
        % Single example case
        dd = cond_time_evolution{2,27};
        n_teor = cond_time_evolution{2,26};
        tt = cond_time_evolution{2,24};
        
        index1 = find(tt == 48);
        
        index = find(cond_hila == num_of_disc);
        index = index+1;
        
        g = cond_time_evolution{index,14};
        n_petrov= cond_time_evolution{index,16};
        error_PGFEM = cond_time_evolution{index,17};
        
        d = cond_time_evolution{index,15};
        n_diff = cond_time_evolution{index,20};
        error_diff = cond_time_evolution{index,21};
        error_diff2 = cond_time_evolution{index,37};
        
        n_petrov_ln = g.*n_petrov;
        n_diff_ln = d.*n_diff;
        n_teor_ln = dd.*n_teor;
        
        semilogx(dd,n_teor_ln(:,1),'k--','LineWidth',2)
        ylim([min(min(n_petrov_ln)),1.05*max(max(n_teor_ln))]);
        hold on
        semilogx(g,n_petrov_ln(:,index1),'b','LineWidth',1.2)
        semilogx(d,n_diff_ln(:,index1),'r','LineWidth',1.2)
        semilogx(dd,n_teor_ln(:,index1),'k--','LineWidth',2)
        semilogx(g,n_petrov_ln(:,end),'b','LineWidth',1.2)
        semilogx(d,n_diff_ln(:,end),'r','LineWidth',1.2)
        semilogx(dd,n_teor_ln(:,end),'k--','LineWidth',2)
        ylabel('Number distribution [d N/d ln(d_p)]')
        lg = legend('Analytical','PGFEM','Sectional','location','best')
        lg.FontSize = 10;
        legend boxoff
        
        
        if jj == 1
            xlabel('Particle diameter [\mum]')
            xlim([1.5e-7,1e-6])
            saveas(gcf,[save_loc,'condensation_test_case_',num2str(length(g)),'.png'])
            h4 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Particle diameter [\mum]'};{'(a)'}])
            xlim([1.5e-7,1e-6])
            subplot(222)
            
        end
        
        semilogy(tt,error_PGFEM,'b','linewidth',1.5)
        hold on
        semilogy(tt,error_diff,'r','linewidth',1.5)
        semilogy(tt,error_diff2,'m','linewidth',1.5)
        ylabel('Relative error [%]')
        lg = legend('PGFEM','Sectional (linear)','Sectional (constant)','location','southeast')
        lg.FontSize = 10;
        legend boxoff
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'condensation_test_case_',num2str(length(g)),'_error.png'])
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(b)'}])
            saveas(gcf,[sl_sub,'fig2'], save_format)
            
            
        end
        
        
    end
    
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Coagulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_case(2) == 1
    close all
    clearvars -except plot_case sl_sub save_format discretization
    load('coagulation_time_evolutions.mat')
    
    if ~exist([pwd,'\figs\figs_coagulation'],'dir')
        mkdir([pwd,'\figs\figs_coagulation'])
    end
    
    save_loc = [pwd,'\figs\figs_coagulation\'];
    
    % Number of elements/bins for the approximation method
    num_of_disc = discretization(2);
    
    % Loop for fetching the variables for the plots
    for ii = 2:size(coagulation_time_evolution_const,1)
        
        % FEM
        
        % number of discretization points
        coag_disc(ii-1,1) = length(coagulation_time_evolution_const{ii,10});
        % average relative error of FEM solution (time step 1 h)
        coag_error_t1(ii-1,1) = coagulation_time_evolution_const{ii,3};
        % Computation times of FEM solution
        coag_time_t1(ii-1,1) = coagulation_time_evolution_const{ii,4};
        % Average relative error of FEM solution (time step 0.1 h)
        coag_error_t2(ii-1,1) = coagulation_time_evolution_const{ii,14};
        % Computation times of FEM solution (0.1 h timesetp)
        coag_time_t2(ii-1,1) = coagulation_time_evolution_const{ii,15};
        
        % Sectional method
        
        % average relative error of sectional method (time step 1 h)
        coag_error_t1_diff(ii-1,1) = coagulation_time_evolution_const{ii,7};
        % Computation times of sectional method (1h timestep)
        coag_time_t1_diff(ii-1,1) = coagulation_time_evolution_const{ii,8};
        % average relative error of sectional solution (time step 0.1 h)
        coag_error_t2_diff(ii-1,1) = coagulation_time_evolution_const{ii,18};
        % Computation times of differenc solution (0.1h timestep)
        coag_time_t2_diff(ii-1,1) = coagulation_time_evolution_const{ii,19};
    end
    
    for jj = 1:2
        
        if jj == 1
            h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            H = fig('width',20,'height',22,'font','Helvetica','fontsize',11)
            subplot(325)
            
        end
        semilogy(coag_disc,coag_error_t1,'bo--','linewidth',1.5)
        hold on
        semilogy(coag_disc,coag_error_t2,'cx--','linewidth',1.5)
        semilogy(coag_disc,coag_error_t1_diff,'rs--','linewidth',1.5)
        semilogy(coag_disc,coag_error_t2_diff,'md--','linewidth',1.5)
        ylabel('Average relative error [%]')
        lg1 = legend('FEM, \Deltat = 1h','FEM, \Deltat = 0.1h','Sec., \Deltat = 1h',...
            'Sec., \Deltat = 0.1h')
        lg1.FontSize = 10;
        legend boxoff
        set(lg1,'Position',[0.262602877695686 ...
            0.238648523944074 0.190476186810032 0.088447650990905])
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        
        %         saveas(gcf,[save_loc,'fig5.fig'])
        
        % This figure is comparison between median error of time evolution to
        % computation time of evolution.
        if jj == 1
            xlabel('Number of bins/elements')
            saveas(gcf,[save_loc,'coagulation_average_error.png'])
            h2 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Number of bins/elements'};{'(e)'}]);
            subplot(326)
            
        end
        loglog(coag_time_t1,coag_error_t1,'bo--','linewidth',1.5)
        hold on
        loglog(coag_time_t2,coag_error_t2,'cx--','linewidth',1.5)
        loglog(coag_time_t1_diff,coag_error_t1_diff,'rs--','linewidth',1.5)
        loglog(coag_time_t2_diff,coag_error_t2_diff,'md--','linewidth',1.5)
        
        ylabel('Average relative error [%]')
        lg2 = legend('FEM, \Deltat = 1h','FEM, \Deltat = 0.1h','Sec., \Deltat = 1h',...
            'Sec., \Deltat = 0.1h')
        lg2.FontSize = 10;
        legend boxoff
        set(lg2,'Position',[0.712339484838544 ...
            0.242256478489529 0.190476186810032 0.088447650990905])
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        xticks([0.01 0.1 1 10 100 1000])
        xticklabels({'0.01','0.1','1','10','100','1000'})
        
        
        % Single test case figures
        index = find(coag_disc == num_of_disc);
        index = index+1;
        g = coagulation_time_evolution_const{index,10};
        n_FEM = coagulation_time_evolution_const{index,12};
        t = coagulation_time_evolution_const{2,20};
        
        d = coagulation_time_evolution_const{index,11};
        n_diff = coagulation_time_evolution_const{index,16};
        
        vv = coagulation_time_evolution_const{2,23};
        n_teor = coagulation_time_evolution_const{2,22};
        % tt = coagulation_time_evolution_const{2,20};
        
        n_FEM_ln = g.*n_FEM;
        n_diff_ln = d.*n_diff;
        n_teor_ln = vv.*n_teor;
        
        index1 = find(t == 48);
        
        if jj == 1
            xlabel('Computation time [s]')
            saveas(gcf,[save_loc,'coagulation_average_error_comp_time.png'])
            h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Computation time [s]'};{'(f)'}])
            subplot(321)
            
        end
        
        semilogx(vv,n_teor_ln(:,1),'k--','linewidth',2)
        hold on
        plot1 = semilogx(g,n_FEM_ln(:,index1),'b','linewidth',1.5);
        plot2 = semilogx(d,n_diff_ln(:,index1),'r','linewidth',1.5);
        plot3 = semilogx(g,n_FEM_ln(:,end),'b','linewidth',1.5);
        plot4 = semilogx(d,n_diff_ln(:,end),'r','linewidth',1.5);
        loglog(vv,n_teor_ln(:,index1),'k--','linewidth',2)
        loglog(vv,n_teor_ln(:,end),'k--','linewidth',2)
        ylim([min(0),1.1*max(n_teor_ln(:,1))])
        xlim([10^(-7),3*10^(-3)])
        ylabel('Number distribution [d N/d ln(v_p)]')
        % plot1.Color(4) = 0.7;
        plot2.Color(4) = 0.6;
        % plot3.Color(4) = 0.7;
        plot4.Color(4) = 0.6;
        % title(['Nodes/bins: ',num2str(length(g)),'. Time step: ', num2str(t(2)-t(1)),'h']);
        lg3 = legend('Analytical','FEM','Sectional','location','northwest')
        lg3.FontSize = 10;
        legend boxoff
        
        error_FEM = coagulation_time_evolution_const{index,13};
        error_diff = coagulation_time_evolution_const{index,17};
        
        if jj == 1
            xlabel('Particle volume [\mum^3]')
            saveas(gcf,[save_loc,'coagulation_test_case_',num2str(length(g)),'.png'])
            h4 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Particle volume [\mum^3]'};{'(a)'}])
            subplot(322)
            
        end
        
        plot(t,error_FEM,'linewidth',1.5)
        hold on
        plot(t,error_diff,'linewidth',1.5)
        
        ylabel('Relative error [%]')
        lg = legend('FEM','Sectional','location','best');
        lg.FontSize = 10;
        legend boxoff
        
        
        %         saveas(gcf,[save_loc,'fig2.fig'])
        
        n_teor_ln = n_teor_ln(vv >= 9e-7,:);
        vv = vv(vv >= 9e-7);
        % Corresponding colormap plots with 100 nodes FEM solution
        for ii = 1:length(t)
            
            n_FEM_inter(:,ii) = interp1(g,n_FEM_ln(:,ii),vv,'linear');
            
        end
        vv1 = vv(vv >= 9e-7);
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'coagulation_test_case_',num2str(length(g)),'_error.png'])
            h5 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11,'border','on');
            labels.xlab = 'Time [h]';
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(b)'}])
            subplot(323)
            labels.xlab = [{'Time [h]'};{'(c)'}];
        end
        % Analytical time evolution
        %         h5 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11,'border','on');
        [d_grid,t_grid] = meshgrid(vv1,t);
        imgrid.d_grid = d_grid;
        imgrid.t_grid = t_grid;
        
        % Labels for the axis
        labels.ylab = 'Particle Volume [\mum^3]';
        labels.yscale = 'Log';
        labels.ytick = [1e-6 1e-5 1e-4 1e-3 1e-2];
        labels.xtick = [0 16 32 48 64 80 96];
        labels.clab = 'Number distribution [d N/d ln(v_p)]';
        labels.size = 13;
        labels.title = [];
        
        % Data
        ImagesForScaling = [n_teor_ln];
        fig1.figno = 5;
        fig1.position = [];
        fig1.subplot = [];
        fig1.clf = 0;
        c = PlotParticleDensityEvolution(n_teor_ln,imgrid,ImagesForScaling,labels,fig1);
        pos2 = get(c,'Position')
        ax = gca;
        ax.FontName = 'Helvetica';
        ax.XLabel.FontSize = 11;
        ax.YLabel.FontSize = 11;
        pos = get(ax,'Position')
        c.Limits = [0,max(max(n_teor_ln))];
        c.Label.FontName = 'Helvetica';
        c.Label.FontSize = 10.5;
        pos2 = [pos2(1)-0.01 pos(2) pos2(3)-0.01 pos(4)];
        c.Position = [pos2];
        pos(3) = pos2(1)-pos(1)-0.005;
        ax.Position = pos;
        ax.XLabel.Position = [get(ax.XLabel,'position')+[17 0 0]]
        
        %         saveas(gcf,[save_loc,'fig3.fig'])
        
        if jj == 1
            saveas(gcf,[save_loc,'coagulation_test_case_colorplot_analytical.png'])
            h6 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11,'border','on');
        elseif jj == 2
            subplot(324)
            labels.xlab = [{'Time [h]'};{'(d)'}];
        end
        
        ImagesForScaling = [n_FEM_inter];
        fig2.figno = 6;
        fig2.position = [];
        fig2.subplot = [];
        fig2.clf = 0;
        c = PlotParticleDensityEvolution(n_FEM_inter,imgrid,ImagesForScaling,labels,fig2);
        pos2 = get(c,'Position')
        ax = gca;
        ax.FontName = 'Helvetica';
        ax.XLabel.FontSize = 11;
        ax.YLabel.FontSize = 11;
        pos = get(ax,'Position')
        c.Limits = [0,max(max(n_teor_ln))];
        c.Label.FontName = 'Helvetica';
        c.Label.FontSize = 10.5;
        pos2 = [pos2(1)-0.01 pos(2) pos2(3)-0.01 pos(4)];
        c.Position = [pos2];
        pos(3) = pos2(1)-pos(1)-0.005;
        ax.Position = pos;
        ax.XLabel.Position = [get(ax.XLabel,'position')+[17 0 0]]
        if jj == 1
            saveas(gcf,[save_loc,'coagulation_test_case_',num2str(length(g)),'_colorplot_FE.png'])
            %     saveas(gcf,[save_loc,'fig4.fig'])
        elseif jj == 2
            saveas(gcf,[sl_sub,'fig3'], save_format)
        end
        
    end
    
    
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_case(3) == 1
    close all
    clearvars -except plot_case sl_sub save_format discretization
    load('GDE_time_evolutions.mat')
    if ~exist([pwd,'\figs\figs_GDE'],'dir')
        mkdir([pwd,'\figs\figs_GDE'])
    end
    save_loc = [pwd,'\figs\figs_GDE\'];
    
    % Number of elements/bins for the approximation method
    num_of_disc = discretization(3);
    
    for ii = 2:size(GDE_time_evolutions,1)
        
        % FE approximation
        
        % Element spacings
        GDE_hila(ii-1,1) = length(GDE_time_evolutions{ii,14});
        
        % FEM average relative errors and computation times
        GDE_error_t1(ii-1,1) = GDE_time_evolutions{ii,3};
        GDE_time_t1(ii-1,1) = GDE_time_evolutions{ii,4};
        GDE_error_t2(ii-1,1) = GDE_time_evolutions{ii,18};
        GDE_time_t2(ii-1,1) = GDE_time_evolutions{ii,19};
        
        % PGFEM average relative errors and computation times
        %     GDE_PGerror_t1(ii-1,1) = GDE_time_evolutions{ii,7};
        %     GDE_PGtime_t1(ii-1,1) = GDE_time_evolutions{ii,8};
        %     GDE_PGerror_t2(ii-1,1) = GDE_time_evolutions{ii,22};
        %     GDE_PGtime_t2(ii-1,1) = GDE_time_evolutions{ii,23};
        
        % Sectional method average relative errors and computation times
        GDE_error_t1_diff(ii-1,1) = GDE_time_evolutions{ii,11};
        GDE_time_t1_diff(ii-1,1) = GDE_time_evolutions{ii,12};
        GDE_error_t2_diff(ii-1,1) = GDE_time_evolutions{ii,26};
        GDE_time_t2_diff(ii-1,1) = GDE_time_evolutions{ii,27};
        
    end
    
    for jj = 1:2 % Standard figures and subplot
        
        if jj == 1
            h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            H = fig('width',20,'height',22','font','Helvetica','fontsize',11)
            subplot(325)
            
        end
        semilogy(GDE_hila,GDE_error_t1,'bo--','linewidth',1.5)
        hold on
        semilogy(GDE_hila,GDE_error_t2,'cx--','linewidth',1.5)
        semilogy(GDE_hila,GDE_error_t1_diff,'rs--','linewidth',1.5)
        semilogy(GDE_hila,GDE_error_t2_diff,'md--','linewidth',1.5)
        yl = ylabel('Average relative error [%]');
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yl.FontSize = 10;
        yticklabels({'0','0.1','1','10','100'})
        lgd = legend('FEM, \Deltat = 1h','FEM, \Deltat = 0.1h','Sec., \Deltat = 1h',...
            'Sec., \Deltat = 0.1h','location','best');
        lgd.FontSize = 9;
        legend boxoff
        
        if jj == 1
            xlabel('Number of bins/elements')
            saveas(gcf,[save_loc,'GDE_average_error.png'])
            h2 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Number of bins/elements'};{'(e)'}]);
            subplot(326)
            
        end
        loglog(GDE_time_t1,GDE_error_t1,'bo--','linewidth',1.5)
        hold on
        loglog(GDE_time_t2,GDE_error_t2,'cx--','linewidth',1.5)
        loglog(GDE_time_t1_diff,GDE_error_t1_diff,'rs--','linewidth',1.5)
        loglog(GDE_time_t2_diff,GDE_error_t2_diff,'md--','linewidth',1.5)
        xlim([min(GDE_time_t1_diff),1200])
        yl = ylabel('Average relative error [%]');
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        yl.FontSize = 10;
        lgd = legend('FEM, \Deltat = 1h','FEM, \Deltat = 0.1h','Sec., \Deltat = 1h',...
            'Sec., \Deltat = 0.1h');
        lgd.FontSize = 8.5;
        legend boxoff
        set(lgd,'Position',[0.731155265930776...
            0.249101629904917 0.173280420126738 0.0788206958024797])
        
        
        if jj == 1
            xlabel('Computation time [s]')
            saveas(gcf,[save_loc,'GDE_average_error_comp_time.png'])
            h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Computation time [s]'};{'(f)'}]);
            subplot(321)
            
        end
        
        % Single example case
        dd = GDE_time_evolutions{2,31};
        n_teor = GDE_time_evolutions{2,30};
        tt = GDE_time_evolutions{2,28};
        
        index1 = find(tt == 48);
        
        index = find(GDE_hila == num_of_disc);
        index = index+1;
        
        g = GDE_time_evolutions{index,14};
        n_FEM= GDE_time_evolutions{index,16};
        error_FEM = GDE_time_evolutions{index,17};
        
        d = GDE_time_evolutions{index,15};
        n_diff = GDE_time_evolutions{index,24};
        error_diff = GDE_time_evolutions{index,25};
        
        n_FEM_ln = g.*n_FEM;
        n_diff_ln = d.*n_diff;
        n_teor_ln = dd.*n_teor;
        
        n_teor(isnan(n_teor)) = 0;
        n_teor_ln(isnan(n_teor_ln)) = 0;
        FEM_particles = trapz(g,n_FEM);
        diff_particles = trapz(d,n_diff);
        teor_particles = trapz(dd,n_teor);
        
        
        FEM_volume = trapz(g,n_FEM_ln);
        diff_volume = trapz(d,n_diff_ln);
        teor_volume = trapz(dd,n_teor_ln);
        
        
        semilogx(dd,n_teor_ln(:,1),'k--','LineWidth',2)
        ylim([min(min(n_FEM_ln)),1.05*max(max(n_teor_ln))]);
        hold on
        semilogx(g,n_FEM_ln(:,index1),'b','LineWidth',1.2)
        semilogx(d,n_diff_ln(:,index1),'r','LineWidth',1.2)
        semilogx(dd,n_teor_ln(:,index1),'k--','LineWidth',2)
        semilogx(g,n_FEM_ln(:,end),'b','LineWidth',1.2)
        semilogx(d,n_diff_ln(:,end),'r','LineWidth',1.2)
        semilogx(dd,n_teor_ln(:,end),'k--','LineWidth',2)
        yl = ylabel('Number distribution [d N/d ln(v_p)]');
        yl.FontSize = 10;
        lg = legend('Analytical','FEM','Sectional','location','best')
        lg.FontSize = 10;
        legend boxoff
        
        
        if jj == 1
            xlabel('Particle volume [\mum^3]')
            %             xlim([1.5e-7,1e-6])
            saveas(gcf,[save_loc,'GDE_test_case_',num2str(length(g)),'.png'])
            h4 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Particle volume [\mum^3]'};{'(a)'}])
            subplot(324)
            
        end
        
        plot(tt,error_FEM,'b','linewidth',1.5)
        hold on
        plot(tt,error_diff,'r','linewidth',1.5)
        yl = ylabel('Relative error [%]');
        yl.FontSize = 10;
        lg = legend('FEM','Sectional','location','northwest')
        lg.FontSize = 10;
        legend boxoff
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'GDE_test_case_',num2str(length(g)),'_error.png'])
            h5 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(d)'}])
            %             saveas(gcf,[sl_sub,'fig2'], save_format)
            subplot(322)
            
        end
        plot(tt,teor_particles,'k--','linewidth',2)
        hold on
        plot(tt,FEM_particles,'b','linewidth',1.2)
        plot(tt,diff_particles,'r','linewidth',1.2)
        lg = legend('Analytical','FEM','Sectional')
        lg.FontSize = 10;
        legend boxoff
        
        yl = ylabel([{'Number concentration of'}; {'particles [cm^{-3}]'}]);
        yl.FontSize = 10;
        
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'GDE_test_case_',num2str(length(g)),'_particle_number.png'])
            h6 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(b)'}])
            subplot(323)
            
        end
        
        plot(tt,teor_volume,'k--','linewidth',2)
        hold on
        plot(tt,FEM_volume,'b','linewidth',1.2)
        plot(tt,diff_volume,'r','linewidth',1.2)
        lg = legend('Analytical','FEM','Sectional','location','best');
        lg.FontSize = 10;
        legend boxoff
        yl = ylabel([{'Volume concentration of'}; {'particles [\mum^3cm^{-3}]'}]);
        yl.FontSize = 10;
        
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'GDE_test_case_',num2str(length(g)),'_particle_volume.png'])
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(c)'}])
            saveas(gcf,[sl_sub,'fig4'], save_format)
            
        end
        
        
        
    end
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Discrete GDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_case(4) == 1
    
    close all
    clearvars -except plot_case sl_sub save_format discretization
    load('Discrete_evolution_GDE_test_case.mat')
    load('discrete_GDE_test_case.mat')
    if ~exist([pwd,'\figs\figs_discrete_GDE'],'dir')
        mkdir([pwd,'\figs\figs_discrete_GDE'])
    end
    save_loc = [pwd,'\figs\figs_discrete_GDE\'];
    
    % Number of elements/bins for the approximation method
    num_of_disc = discretization(4);
    
    for ii = 2:size(discrete_GDE_evolutions_numerical,1)
        
        % FE approximation
        
        % Element spacings
        GDE_hila(ii-1,1) = length(discrete_GDE_evolutions_numerical{ii,14});
        
        % FEM average relative errors and computation times
        GDE_error_t1(ii-1,1) = discrete_GDE_evolutions_numerical{ii,3};
        GDE_time_t1(ii-1,1) = discrete_GDE_evolutions_numerical{ii,4};
        
        
        % PGFEM average relative errors and computation times
        GDE_PGerror_t1(ii-1,1) = discrete_GDE_evolutions_numerical{ii,7};
        GDE_PGtime_t1(ii-1,1) = discrete_GDE_evolutions_numerical{ii,8};
        
        
        % Sectional method average relative errors and computation times
        GDE_error_t1_diff(ii-1,1) = discrete_GDE_evolutions_numerical{ii,11};
        GDE_time_t1_diff(ii-1,1) = discrete_GDE_evolutions_numerical{ii,12};
        
        
    end
    
    for jj = 1:2 % Standard figures and subplot
        
        if jj == 1
            h1 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            H = fig('width',20,'height',22','font','Helvetica','fontsize',11)
            subplot(325)
            
        end
        loglog(GDE_hila,GDE_error_t1,'bo-','linewidth',1.5)
        hold on
        loglog(GDE_hila,GDE_PGerror_t1,'ro-','linewidth',1.5)
        loglog(GDE_hila,GDE_error_t1_diff,'co-','linewidth',1.5)
        yl = ylabel('Average relative error [%]');
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yl.FontSize = 10;
        yticklabels({'0','0.1','1','10','100'})
        lgd = legend('FEM','PGFEM','Sectional','location','best');
        lgd.FontSize = 9;
        legend boxoff
        
        if jj == 1
            xlabel('Number of bins/elements')
            saveas(gcf,[save_loc,'discrete_GDE_average_error.png'])
            h2 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Number of bins/elements'};{'(e)'}]);
            subplot(326)
            
        end
        loglog(GDE_time_t1,GDE_error_t1,'bo-','linewidth',1.5)
        hold on
        loglog(GDE_PGtime_t1,GDE_PGerror_t1,'ro-','linewidth',1.5)
        loglog(GDE_time_t1_diff,GDE_error_t1_diff,'co-','linewidth',1.5)
        
        yl = ylabel('Average relative error [%]');
        ylim([0,100])
        yticks([0 0.1 1 10 100])
        yticklabels({'0','0.1','1','10','100'})
        yl.FontSize = 10;
        lgd = legend('FEM','PGFEM','Sectional','location','best');
        lgd.FontSize = 9;
        legend boxoff
        
        
        if jj == 1
            xlabel('Computation time [s]')
            saveas(gcf,[save_loc,'discrete_GDE_average_error_comp_time.png'])
            h3 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
        elseif jj == 2
            xlabel([{'Computation time [s]'};{'(f)'}]);
            subplot(321)
            
        end
        
        
        index = find(GDE_hila == num_of_disc);
        index = index+1;
        
        g = discrete_GDE_evolutions_numerical{index,14};
        n_FEM= discrete_GDE_evolutions_numerical{index,1};
        error_FEM = discrete_GDE_evolutions_numerical{index,2};
        n_PGFEM= discrete_GDE_evolutions_numerical{index,5};
        error_PGFEM = discrete_GDE_evolutions_numerical{index,6};
        
        d = discrete_GDE_evolutions_numerical{index,15};
        n_diff = discrete_GDE_evolutions_numerical{index,9};
        error_diff = discrete_GDE_evolutions_numerical{index,10};
        
        n_FEM_ln = g.*n_FEM;
        n_PGFEM_ln = g.*n_PGFEM;
        n_diff_ln = d.*n_diff;
        n_disc_ln = d_disc.*n_disc;
        
        n_disc(isnan(n_disc)) = 0;
        n_disc_ln(isnan(n_disc_ln)) = 0;
        FEM_particles = trapz(g,n_FEM);
        PGFEM_particles = trapz(g,n_PGFEM);
        diff_particles = trapz(d,n_diff);
        disc_particles = trapz(d_disc,n_disc);
        
        
        FEM_volume = trapz(g,n_FEM_ln);
        PGFEM_volume = trapz(g,n_PGFEM_ln);
        diff_volume = trapz(d,n_diff_ln);
        disc_volume = trapz(d_disc,n_disc_ln);
        
        
        semilogx(d_disc,n_disc_ln(:,1),'k--','LineWidth',2)
        ylim([min(min(n_FEM_ln)),1.05*max(max(n_disc_ln))]);
        xlim([10^(-9),10^(-5)])
        xticks([1e-9 1e-8 1e-7 1e-6 1e-5])
        hold on
        semilogx(g,n_FEM_ln(:,end),'b','LineWidth',1.2)
        semilogx(g,n_PGFEM_ln(:,end),'r','LineWidth',1.2)
        semilogx(d,n_diff_ln(:,end),'c','LineWidth',1.2)
        semilogx(d_disc,n_disc_ln(:,end),'k--','LineWidth',2)
        yl = ylabel('Number distribution [d N/d ln(v_p)]');
        yl.FontSize = 10;
        lg = legend('Discrete','FEM','PGFEM','Sectional','location','best')
        lg.FontSize = 9;
        legend boxoff
        
        
        if jj == 1
            xlabel('Particle volume [\mum^3]')
            
            %             xlim([1.5e-7,1e-6])
            saveas(gcf,[save_loc,'discrete_GDE_test_case_',num2str(length(g)),'.png'])
            h4 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Particle volume [\mum^3]'};{'(a)'}])
            subplot(324)
            
        end
        
        tt = tt_disc;
        
        plot(tt,error_FEM,'b','linewidth',1.5)
        hold on
        plot(tt,error_PGFEM,'r','linewidth',1.5)
        plot(tt,error_diff,'c','linewidth',1.5)
        yl = ylabel('Relative error [%]');
        yl.FontSize = 10;
        lg = legend('FEM','PGFEM','Sectional','location','northwest')
        lg.FontSize = 9;
        legend boxoff
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'discrete_GDE_test_case_',num2str(length(g)),'_error.png'])
            h5 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(d)'}])
            %             saveas(gcf,[sl_sub,'fig2'], save_format)
            subplot(322)
            
        end
        plot(tt,disc_particles,'k--','linewidth',2)
        hold on
        plot(tt,FEM_particles,'b','linewidth',1.2)
        plot(tt,PGFEM_particles,'r','linewidth',1.2)
        plot(tt,diff_particles,'c','linewidth',1.2)
        lg = legend('Discrete','FEM','PGFEM','Sectional')
        lg.FontSize = 9;
        legend boxoff
        if jj == 2
            set(lg,'Position',[0.783952768725903 0.83773851102879 0.120370368399317 0.0835842780091546])
        end
        yl = ylabel([{'Number concentration of'}; {'particles [cm^{-3}]'}]);
        yl.FontSize = 10;
        
        
        
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'discrete_GDE_test_case_',num2str(length(g)),'_particle_number.png'])
            h6 = fig('width',12.75,'height',12,'font','Helvetica','fontsize',11);
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(b)'}])
            subplot(323)
            
        end
        
        plot(tt,disc_volume,'k--','linewidth',2)
        hold on
        plot(tt,FEM_volume,'b','linewidth',1.2)
        plot(tt,FEM_volume,'r','linewidth',1.2)
        plot(tt,diff_volume,'c','linewidth',1.2)
        lg = legend('Discrete','FEM','PGFEM','Sectional','location','best');
        lg.FontSize = 9;
        legend boxoff
        yl = ylabel([{'Volume concentration of'}; {'particles [\mum^3cm^{-3}]'}]);
        yl.FontSize = 10;
        
        if jj == 1
            xlabel('Evolution time [h]')
            saveas(gcf,[save_loc,'discrete_GDE_test_case_',num2str(length(g)),'_particle_volume.png'])
            
            
        elseif jj == 2
            xlabel([{'Evolution time [h]'};{'(c)'}])
            saveas(gcf,[sl_sub,'fig5'], save_format)
            
        end
        
        
        
    end
    
    d = [min(d_disc);d;max(d_disc)];
    n_diff_ln = [n_diff_ln(1,:);n_diff_ln;n_diff_ln(end,:)];
    for ii = 1:length(tt_disc)
        
        n_PGFEM_inter(:,ii) = interp1(g,n_PGFEM_ln(:,ii),d_disc,'linear');
        n_diff_inter(:,ii) = interp1(d,n_diff_ln(:,ii),d_disc,'linear');
        
    end
    
    
    H = fig('width',20,'height',13','font','Helvetica','fontsize',11)
    h1 = subplot(131);
    [d_grid,t_grid] = meshgrid(d_disc,tt_disc);
    imgrid.d_grid = d_grid;
    imgrid.t_grid = t_grid;
    % labels
    labels.xlab = [{'Time [s]'};{'(a)'}];
    labels.ylab = 'Particle Volume [\mum^3]';
    labels.yscale = 'Log';
    labels.ytick = [1e-8 1e-7 1e-6 1e-5];
    % labels.xtick = [0 10 20 30 40 48];
    labels.xtick = [0 0.25 0.5];
    labels.clab = 'Number distribution [d N/d ln(v_p)]';
    labels.size = 10.5;
    labels.title = [];
    % Data
    ImagesForScaling = [n_disc_ln];
    fig1.figno = 7;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;
    c = PlotParticleDensityEvolution(n_disc_ln,imgrid,ImagesForScaling,labels,fig1);
    pos2 = get(c,'Position')
    ax = gca;
    ax.FontName = 'Helvetica';
    pos1 = ax.Position;
    ax.Position = [pos1(1),pos1(2)+0.05,pos1(3),pos1(4)-0.05];
    pos = get(ax,'Position')
    c.Limits = [0,max(max(n_disc_ln))];
    c.Label.FontName = 'Helvetica';
    c.Label.FontSize = 10;
    pos2 = [pos2(1)-0.01 pos(2) pos2(3)-0.005 pos(4)];
    c.Position = [pos2];
    pos(3) = pos2(1)-pos(1)-0.005;
    ax.Position = pos;
    % ax.XLabel.Position = [get(ax.XLabel,'position')+[17 0 0]]
    
    h2 = subplot(132);
    [d_grid,t_grid] = meshgrid(d_disc,tt_disc);
    imgrid.d_grid = d_grid;
    imgrid.t_grid = t_grid;
    % labels
    labels.xlab = [{'Time [s]'};{'(b)'}];
    labels.ylab = 'Particle Volume [\mum^3]';
    labels.yscale = 'Log';
    labels.ytick = [1e-8 1e-7 1e-6 1e-5];
    % labels.xtick = [0 10 20 30 40 48];
    labels.xtick = [0 0.25 0.5];
    labels.clab = 'Number distribution [d N/d ln(v_p)]';
    labels.size = 10.5;
    labels.title = [];
    % Data
    ImagesForScaling = [n_PGFEM_inter];
    fig1.figno = 7;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;
    c = PlotParticleDensityEvolution(n_PGFEM_inter,imgrid,ImagesForScaling,labels,fig1);
    pos2 = get(c,'Position')
    ax = gca;
    ax.FontName = 'Helvetica';
    pos1 = ax.Position;
    ax.Position = [pos1(1),pos1(2)+0.05,pos1(3),pos1(4)-0.05];
    pos = get(ax,'Position')
    c.Limits = [0,max(max(n_disc_ln))];
    c.Label.FontName = 'Helvetica';
    c.Label.FontSize = 10;
    pos2 = [pos2(1)-0.01 pos(2) pos2(3)-0.005 pos(4)];
    c.Position = [pos2];
    pos(3) = pos2(1)-pos(1)-0.005;
    ax.Position = pos;
    
    % ax.XLabel.Position = [get(ax.XLabel,'position')+[17 0 0]]
    h3 = subplot(133);
    [d_grid,t_grid] = meshgrid(d_disc,tt_disc);
    imgrid.d_grid = d_grid;
    imgrid.t_grid = t_grid;
    % labels
    labels.xlab = [{'Time [s]'};{'(c)'}];
    labels.ylab = 'Particle Volume [\mum^3]';
    labels.yscale = 'Log';
    labels.ytick = [1e-8 1e-7 1e-6 1e-5];
    % labels.xtick = [0 10 20 30 40 48];
    labels.xtick = [0 0.25 0.5];
    labels.clab = 'Number distribution [d N/d ln(v_p)]';
    labels.size = 10.5;
    labels.title = [];
    % Data
    ImagesForScaling = [n_diff_inter];
    fig1.figno = 7;
    fig1.position = [];
    fig1.subplot = [];
    fig1.clf = 0;
    c = PlotParticleDensityEvolution(n_diff_inter,imgrid,ImagesForScaling,labels,fig1);
    pos2 = get(c,'Position')
    ax = gca;
    ax.FontName = 'Helvetica';
    pos1 = ax.Position;
    ax.Position = [pos1(1),pos1(2)+0.05,pos1(3),pos1(4)-0.05];
    pos = get(ax,'Position')
    c.Limits = [0,max(max(n_disc_ln))];
    c.Label.FontName = 'Helvetica';
    c.Label.FontSize = 10;
    pos2 = [pos2(1)-0.01 pos(2) pos2(3)-0.005 pos(4)];
    c.Position = [pos2];
    pos(3) = pos2(1)-pos(1)-0.005;
    ax.Position = pos;
    saveas(gcf,[sl_sub,'fig6'], save_format)
    % ax.XLabel.Position = [get(ax.XLabel,'position')+[17 0 0]]
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
