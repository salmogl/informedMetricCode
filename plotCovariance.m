function plotCovariance(covariance,title)

sizeRows    = length(covariance);
fontS       = 20;
coorZoom    = [750 900];
as          = -46; % View parameters
el          = 11; 


figure;
mesh(1:sizeRows,1:sizeRows,covariance);

zlim([-2 30])
xlim([200 sizeRows])
ylim([200 sizeRows])

% xlim([1 sizeRows]);ylim([1 sizeRows])

xlabel('$i$','FontSize',fontS)
ylabel('$j$','FontSize',fontS)
zlabel(title,'FontSize',fontS)
% create a new pair of axes inside current figure
axes('position',[.13 .5 .3 .35])
box on % put box around new pair of axes
mesh(1:sizeRows,1:sizeRows,covariance) % plot on new axes
% view(2);
% axis tight


set(gcf,'Renderer','Zbuffer')
xlim(coorZoom)
ylim(coorZoom)
zlim([11 28])
view([as el])

% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% set(gca,'ztick',[])
% set(gca,'zticklabel',[])
set(gca,'color','none')
axis off

