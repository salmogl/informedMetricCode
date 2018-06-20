
load('C:\Users\salmogl\Documents\MATLAB\GenesExspressrion\results\workspace_gapStatProjNorm7Clust');
knn_10_10_60    = knn;
sumD_10_10_60   = sumD;

load('C:\Users\salmogl\Documents\MATLAB\GenesExspressrion\results\workspace_gapStatProjNorm7Clust_10_2_60');
knn_10_2_60     = knn;
sumD_10_2_60    = sumD;

knn             = [knn_10_10_60 knn_10_2_60];
[knn,indSort]   = sort(knn);
sumD            = [sumD_10_10_60 sumD_10_2_60];
sumD            = sumD(indSort);

%%
figure;
plot(knn,sumD,'k','linewidth',2)
hold on
% plot(knn(1:end-1),sumD(2:end)+sd_k(2:end))
% hold on
plot(knn(14),sumD(14),'o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor'...
    ,0.7*[1 1 1],'MarkerSize',8);
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('Gap Statistic, $g(N,7)$','Interpreter','Latex','FontSize',fSize)
xlim([10 60])
text(30,3.6,'$g(N=34,7)$','Interpreter','latex','FontSize',fSize-4)