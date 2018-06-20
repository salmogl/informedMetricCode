function plot3dData(data,color,title_)

pointSize = 75;
scatter3(data(1,:),data(2,:),data(3,:),pointSize,color, 'Fill'); 

xlabel('$\varphi_1$','Interpreter','latex','FontSize',15)
ylabel('$\varphi_2$','Interpreter','latex','FontSize',15)
zlabel('$\varphi_3$','Interpreter','latex','FontSize',15)
h=title(title_);
set(h,'Interpreter','latex','FontSize',15)
