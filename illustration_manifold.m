clear all;
close all;
clc;

[x,y] = meshgrid(0:.05:8,0:.05:8);

X = ((1/2)*(x-(1/20)*(y-4).^2)+1).*cos(x);
Y = y;
Z = ((1/2)*(x-(1/20)*(y-4).^2)+1).*sin(x);

tanPointx=[4,5,6.8,7.5];
tanPointy=[7.5,0.5,4,4.5];


dXdx = (1/2).*cos(tanPointx)-((1/2)*(tanPointx-(1/20)*(tanPointy-4).^2)+1).*sin(tanPointx);
dYdx = zeros(size(tanPointx));
dZdx = (1/2).*sin(tanPointx)+((1/2)*(tanPointx-(1/20)*(tanPointy-4).^2)+1).*cos(tanPointx);

dXdy = -(1/2)*(2*(1/20)*(tanPointy-4)).*cos(tanPointx);
dYdy = ones(size(tanPointx));
dZdy = -(1/2)*(2*(1/20)*(tanPointy-4)).*sin(tanPointx);

basisx=[dXdx;dYdx;dZdx];
basisx=basisx*diag(1./sqrt(diag(basisx'*basisx)));
basisy=[dXdy;dYdy;dZdy];
basisy=basisy*diag(1./sqrt(diag(basisy'*basisy)));


X = ((1/2)*(x-(1/20)*(y-4).^2)+1).*cos(x);
Y = y;
Z = ((1/2)*(x-(1/20)*(y-4).^2)+1).*sin(x);

tanPointX = ((1/2)*(tanPointx-(1/20)*(tanPointy-4).^2)+1+0.01).*cos(tanPointx);
tanPointY = tanPointy;
tanPointZ = ((1/2)*(tanPointx-(1/20)*(tanPointy-4).^2)+1+0.01).*sin(tanPointx);

plane1=repmat([tanPointX(1);tanPointY(1);tanPointZ(1)],[1,4])+2*[-basisx(:,1)-basisy(:,1),-basisx(:,1)+basisy(:,1),+basisx(:,1)+basisy(:,1),+basisx(:,1)-basisy(:,1)];
plane2=repmat([tanPointX(2);tanPointY(2);tanPointZ(2)],[1,4])+2*[-basisx(:,2)-basisy(:,2),-basisx(:,2)+basisy(:,2),+basisx(:,2)+basisy(:,2),+basisx(:,2)-basisy(:,2)];
plane3=repmat([tanPointX(3);tanPointY(3);tanPointZ(3)],[1,4])+2*[-basisx(:,3)-basisy(:,3),-basisx(:,3)+basisy(:,3),+basisx(:,3)+basisy(:,3),+basisx(:,3)-basisy(:,3)];

tri=delaunay(x,y);
tri=tri'; % reshape in 3 x ntri
patch(X(tri),Y(tri),Z(tri),x(tri),'FaceAlpha',0.7);

axis('equal');
shading interp;
hold on;
scatter3(tanPointX,tanPointY,tanPointZ,16,'MarkerEdgeColor','k','MarkerFaceColor','r');

view(-45,16); %# sets the perspective

axis([-8,8,-3,13,-8,8]/1.3);

patch(plane1(1,:),plane1(2,:),plane1(3,:),'r','FaceAlpha',0.7);
patch(plane2(1,:),plane2(2,:),plane2(3,:),'r','FaceAlpha',0.7);
patch(plane3(1,:),plane3(2,:),plane3(3,:),'r','FaceAlpha',0.7);
% annotation(figure1,'textbox',...
%     [0.374214285714286 0.42857142857143 0.0615 0.061904761904762],'String','y',...
%     'FitBoxToText','off');
% str = '$$y_1$$';
% text(tanPointX(1),tanPointY(1),tanPointZ(1),str,'Interpreter','latex')
% str = '$$y_2$$';
% text(tanPointX(2),tanPointY(2),tanPointZ(2),str,'Interpreter','latex')
% str = '$$y_3$$';
% text(tanPointX(3),tanPointY(3),tanPointZ(3),str,'Interpreter','latex')
factor=1.2;
annotation('textbox',[factor*(0.43-0.5)+0.5,factor*(0.37-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$y_1$$'},'interpreter','latex','LineStyle','none');
annotation('textbox',[factor*(0.57-0.5)+0.5,factor*(0.315-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$y_2$$'},'interpreter','latex','LineStyle','none');
annotation('textbox',[factor*(0.556-0.5)+0.5,factor*(0.553-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$y_3$$'},'interpreter','latex','LineStyle','none');

annotation('textbox',[factor*(0.393-0.5)+0.5,factor*(0.39-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$T_{\mathbf{y_1}}\mathcal{M}$$'},'interpreter','latex','LineStyle','none');
annotation('textbox',[factor*(0.595-0.5)+0.5,factor*(0.285-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$T_{\mathbf{y_2}}\mathcal{M}$$'},'interpreter','latex','LineStyle','none');
annotation('textbox',[factor*(0.558-0.5)+0.5,factor*(0.632-0.5)+0.5,.1,.1],'FontSize',12,'String',{'$$T_{\mathbf{y_3}}\mathcal{M}$$'},'interpreter','latex','LineStyle','none');

