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

xlabel('$i$','FontSize',fontS)
ylabel('$j$','FontSize',fontS)
zlabel(title,'FontSize',fontS)
axes('position',[.13 .5 .3 .35])
box on
mesh(1:sizeRows,1:sizeRows,covariance)

set(gcf,'Renderer','Zbuffer')
xlim(coorZoom)
ylim(coorZoom)
zlim([11 28])
view([as el])

set(gca,'color','none')
axis off

