load('workspace_k10_p5_knn_10_5_40');
pClustAllK_10_5_40 = pClustAllK(1:7);
pPcaAllK_10_5_40   = pPcaAllK(1:7);
pVarClust_10_5_40  = pVarClust(1:7);

load('workspace_k10_p5_knn_45_5_82');
pClustAllK_45_5_82 = pClustAllK;
pPcaAllK_45_5_82   = pPcaAllK;
pVarClust_45_5_82  = pVarClust;


pClustAllK = [pClustAllK_10_5_40 pClustAllK_45_5_82];
pVarClust  = [pVarClust_10_5_40  pVarClust_45_5_82];
pPcaAllK   = [pPcaAllK_10_5_40  pPcaAllK_45_5_82];

knn = [10:5:75 82];

fSize = 18;
std0 = sqrt(pVarClust);
figure;
% plot(x,y,'r.-');
h(1) = plot(knn,pPcaAllK,'r.-')
hold on
h(2) = patch([knn fliplr(knn)],[pClustAllK+std0 fliplr(pClustAllK-std0)],0.8*[1 1 1]);
h(3) = plot(knn,pClustAllK,'b.-')


% hold on
% plot(knn,pClustAllK + sqrt(pVarClust),'--b')
% plot(knn,pClustAllK - sqrt(pVarClust),'--b')
% le = legend('Clustering ($STD$)','Clustering','PCA')
le = legend(h([1 3]),'$p_{LM}$','$p_{ILM}$')

set(le,'Interpreter','Latex','FontSize',fSize-3)
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('P-value','Interpreter','Latex','FontSize',fSize)
xlim([10 65])
ylim([0 0.5])