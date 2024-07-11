function [c] = PlotParticleDensityEvolution(f_image,imgrid,ImagesForScaling,labels,fig,position);
% Colorplot plot with the logscale y axis

% if isempty(imgrid)
%     [d_grid,t_grid] = meshgrid(d,t);
%     imgrid.d_grid = d_grid;
%     imgrid.t_grid = t_grid;
% end

% Matthew Ozon
% University of Eastern Finland
% Department of Applied Physics
 

minval = min(min(ImagesForScaling));
maxval = max(max(ImagesForScaling));

% figure(fig.figno), 

if fig.clf
  clf
end
if ~isempty(fig.position)
  set(fig.figno,'Position',fig.position)
end
if ~isempty(fig.subplot)
  subplot(fig.subplot(1),fig.subplot(2),fig.subplot(3))
end    

pl = pcolor(imgrid.t_grid,imgrid.d_grid,f_image');
set(pl,'EdgeColor','None')
set(gca,'YScale',labels.yscale)
set(gca,'YTick',labels.ytick)
c = colorbar;
c.Label.String = labels.clab;
c.Label.FontSize = labels.size;

xlabel(labels.xlab,'fontsize',labels.size), ylabel(labels.ylab,'fontsize',labels.size),
set(gca,'CLim',[minval,maxval])
set(gca,'XTick',labels.xtick)
title(labels.title)

