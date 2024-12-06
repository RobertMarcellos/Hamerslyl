%% This function is used to draw the test functions
clearvars -except sl_sub save_format 
close all
clc

% Node points for single element
nodes = [-0.5,0,0.5]';

xx = linspace(-1,1,1000)';
xx = sort([xx;0]);

% Calling function from below
phi = Galerkin(nodes,xx);

h1 = fig('width',20,'height',10,'font','Helvetica','fontsize',11);
subplot(121)
plot(xx,phi,'linewidth',1.5)
ylim([-1,2])
xticks([-0.5 0 0.5])
xticklabels({'v_{i-1}','v_{i}','v_{i+1}'})
yticks([-1 0 1 2])
grid on

subplot(122)
plot(xx,Petrov(nodes,xx,0.3),'linewidth',1.5)
hold on
plot(xx,Petrov(nodes,xx,0.5),'linewidth',1.5)
plot(xx,Petrov(nodes,xx,0.7),'linewidth',1.5)
legend('\epsilon = 0.3','\epsilon = 0.5','\epsilon = 0.7')
legend boxoff
xticks([-0.5 0 0.5])
xticklabels({'v_{i-1}','v_{i}','v_{i+1}'})
yticks([-1 0 1 2])
grid on

saveas(gcf,[sl_sub,'fig1'], save_format)






% Function for the standard Galerkin basis/test function
% Equation (7)
% Standard Galerkin basis/test function
function yy = Galerkin(nodes,xx)

h = nodes(2)-nodes(1);

for i = 1:length(xx)
    
    if xx(i) >= nodes(1) && xx(i) <= nodes(2)
        
        yy(i,1) = (1/h)*(xx(i) - nodes(1)); 
        
    elseif xx(i) > nodes(2) && xx(i) <= nodes(3)
        yy(i,1) = (1/h)*(-xx(i) + nodes(3)); 
        
    else
        yy(i,1) = 0;
        
    end

end
end

% Petrov-Galerkin test function (Equations 11 and 12)
function yy = Petrov(nodes,xx,UF)
h = nodes(2)-nodes(1);
for i = 1:length(xx)

    if xx(i) >= nodes(1) && xx(i) <= nodes(2)
        
        yy(i,1) = Galerkin(nodes,xx(i)) + (3/2)*UF*(4/h^2)*(xx(i)-nodes(1))*(nodes(2)-xx(i));
        
    elseif xx(i) > nodes(2) && xx(i) <= nodes(3)
        yy(i,1) = Galerkin(nodes,xx(i)) -(3/2)*UF*(4/h^2)*(xx(i)-nodes(2))*(nodes(3)-xx(i));
        
    else
        yy(i,1) = 0;
        
    end

end
end