clear all;
close all hidden;
clc; 

addpath('..\..\localMahalanobis')
addpath('..\..\mahalanobisEquivalence\3DQuest\3DQuest\Questionnaire')
addpath('kmplot')
addpath('logrank')
addpath('..\')

load('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_CANDF.mat');

plotFlag = false;

if (plotFlag)
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

pSPARCoC = 0.0096; % p-Value obtained by SPARCoC methode
data = dataStruct.genesData;
month = dataStruct.month;
vitalStat = dataStruct.vitalStat;

[val,indMaxStd] = sort(std(data,[],2),'descend');
data = data(indMaxStd(2:201),:);
[M,N] = size(data);
 

rank = 6;
nClust = 7;
% [class,centroids] = kmeans(data,nClust,'MaxIter',1e6);
[class,centroids] = kmeans(data,nClust,'MaxIter',1e6,'OnlinePhase','on');
[val,indSort] = sort(class);

figure;
% ratio = 8;
fSize = 15;
% subplot(1,ratio,1:ratio-1)
HMObject = HeatMap(data(indSort,:));
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)

% colormap(1-gray)

% subplot(1,ratio,ratio)
figure;
imagesc(class(indSort,:))
title('Clusters','Fontsize',fSize)
set(gca,'Xtick',[],'Ytick',[])
colormap gray


figure;
% subplot(1,ratio,1:ratio-1)
HMObject = HeatMap(data);
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)


figure;
imagesc(cov(data(indSort,:)'))
title('Subjects Covariance')

%%
kIter       = 1; % PC of people is avaraged along these iterations
pIter       = 1; % p-value is avaraged along these iterations

knn         = 30; %[45:5:75 82];%[10 15 20 25 30 35 40 45]; % 20 is the optimal
pClustAllP  = zeros(1,pIter);
pPcaAllP    = zeros(1,pIter);
pClustAllK  = zeros(1,length(knn));
pVarClust   = zeros(1,length(knn));
pPcaAllK    = zeros(1,length(knn));
sumD        = zeros(1,length(knn));
sd_k        = zeros(1,length(knn));

args.knn        = knn;
args.replicates = 300;
args.nClust     = nClust;
args.rank       = rank;


maxMonth        = 60;

for mm = 1:length(knn)
    knn(mm)
    args.knn    = knn(mm);

for pp = 1:pIter
    
    votesClust  = zeros(N,kIter);
    VclustAll   = zeros(N,kIter);
    
parfor kk = 1:kIter
    
[mDistPca,mDistClust] = PCA_Reconstruction(data,args);
% mDistPca = squareform(pdist(data'));
% mDistClust = squareform(pdist(data'));
% load('mDistances.mat')

%%
eps = 1;%2
[~, nn_dist] = knnsearch(data', data', 'k', 20);
sigma = eps * median(mDistPca(:));

affPca = exp(-mDistPca.^2/(2*sigma^2)); 
affPca = stochastic(affPca);
eigsnum = 80;
[eigvecs, eigvals] = eigs(affPca, eigsnum);
embeddingPca= eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);


sigma = eps * median(mDistClust(:));
affClust = exp(-mDistClust.^2/(2*sigma^2)); 
affClust = stochastic(affClust);
[eigvecs, eigvals] = eigs(affClust, eigsnum);
embeddingClust = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);



color = jet(max(month)+1);
colorsVital = color(month,:);
colorsVital(vitalStat == 1,:) = repmat(color(max(month)+1,:),length(find(vitalStat)), 1);
colorsStage = dataStruct.stage; %mod(dataStruct.stage,4); %dataStruct.stage; %floor(dataStruct.stage/4);
% [val,deadOrder] = sort(month(~vitalStat));
% colorsVital0 = jet(max(month(~vitalStat)));
% colorsVital0 = colorsVital0(month(~vitalStat),:);



% colorsVital0 = jet(length(month(~vitalStat)));
% colorsVital0 = colorsVital0(deadOrder,:);

% colorsVital0 = jet(round((max(month(~vitalStat)))));
% colorsVital0 = colorsVital0(round((month(~vitalStat))),:);
% deadOrder    = 1:length(find(~vitalStat));

[val,deadOrder] = sort(month(~vitalStat));
colorsVital0 = jet(length(month(~vitalStat)));
[~,kLargest] = sort(deadOrder);
colorsVital0 = colorsVital0(kLargest,:);

deadIdx = 1:N;
deadIdx = deadIdx(~vitalStat);
[~,~,Vclust]    = svd(embeddingClust');


if (plotFlag)
    figure;
    subplot(221)
    plot3dData(embeddingPca(:,1:3)',colorsStage,'Embedding Colored By Stage (PCA)')
    subplot(222)
    plot3dData(embeddingPca(:,1:3)',colorsVital,'Embedding Colored By Month (PCA)')
    subplot(223)
    plot3dData(embeddingPca(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (PCA)')
    colorbar;
    subplot(224)
    plot3dData(embeddingClust(:,1:3)',vitalStat,'Embedding Colored By Vital Status (Dead) (PCA)')
    colorbar;

    figure;
    subplot(221)
    plot3dData(embeddingClust(:,1:3)',colorsStage,'Embedding Colored By Stage (Clustering)')
    subplot(222)
    plot3dData(embeddingClust(:,1:3)',colorsVital,'Embedding Colored By Month (Clustering)')
    subplot(223)
    plot3dData(embeddingClust(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (Clustering)')
    colorbar;
    subplot(224)
    plot3dData(embeddingClust(:,1:3)',vitalStat,'Embedding Colored By Vital Status (Dead) (Clustering)')
    colorbar;
    
    
    figure;
    subplot(221)
    plot3dData(Vclust(:,1:3)',colorsStage,'Embedding Colored By Stage (Clustering)')
    subplot(222)
    plot3dData(Vclust(:,1:3)',colorsVital,'Embedding Colored By Month (Clustering)')
    subplot(223)
    plot3dData(Vclust(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (Clustering)')
    colorbar;
    subplot(224)
    plot3dData(Vclust(:,1:3)',vitalStat,'Embedding Colored By Vital Status (Dead) (Clustering)')
    colorbar;
    
    figure;
    plot3dData(embeddingClust(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (Clustering)')

    figure;
    plot3dData(Vclust(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (Clustering)')
    
    [theta,rho,z] = cart2pol(embeddingClust(:,1),embeddingClust(:,2),embeddingClust(:,3));
    polarEmbedding = [theta, rho, z];
    figure;
    plot3dData(polarEmbedding(deadIdx,1:3)',colorsVital0,'Embedding Colored By Month (Dead) (Clustering)')

    

end
% subplot(223)
% plot3dData(dataProjClust(1:3,:),colorsStage,'First 3 PC Colored By Stage (Clustering)')
% subplot(224)
% plot3dData(dataProjClust(1:3,:),colorsVital,'First 3 PC Colored By Month (Clustering)')

%%



% groupsPca   = kmeans(embeddingPca,2,'MaxIter',1e6,'OnlinePhase','on','Replicates',50);
% groupsClust   = kmeans(embeddingClust,2,'MaxIter',1e6,'OnlinePhase','on','Replicates',50);

% groupsPca = kmeans(embeddingPca,2,'MaxIter',1e6);
% groupsClust = kmeans(embeddingClust,2,'MaxIter',1e6);

% groupsPca = KmeansSvd(embeddingPca',2,true,plotFlag,'PCA');
% groupsClust = KmeansSvd(embeddingClust',2,true,plotFlag,'Clustering');

[~,~,Vclust]    = svd(embeddingClust');
% trytry
% Vclust          = embeddingClust(:,1);
%
groupsClust     = ones(N,1);
Vclust          = Vclust(:,1)*sign(Vclust(1,1));
groupsClust(Vclust > 0) = 2;

[~,~,Vpca]      = svd(embeddingPca');
% trytry
% Vpca            = embeddingPca(:,1);
%
groupsPca       = ones(N,1);
Vpca            = Vpca(:,1)*Vpca(1,1);
groupsPca(Vpca > 0) = 2;


[vHistPca0, tAxisPca0] = getVitalHist(vitalStat(groupsPca == 1),month(groupsPca == 1));
[vHistPca1, tAxisPca1] = getVitalHist(vitalStat(groupsPca == 2),month(groupsPca == 2));
% [vHistPca2, tAxisPca2] = getVitalHist(vitalStat(groupsPca == 3),month(groupsPca == 3));

% figure;
% [~,~,a1,b1] = kmplot([ month(groupsPca == 1) , vitalStat(groupsPca == 1) ]);
% [~,~,a2,b2] = kmplot([ month(groupsPca == 2) , vitalStat(groupsPca == 2) ]);
% [~,~,a3,b3] = kmplot([ month(groupsPca == 3) , vitalStat(groupsPca == 3) ]);

if (plotFlag)
    figure;
    subplot(121)
    stairs(tAxisPca0,vHistPca0,'r')
    hold on
    stairs(tAxisPca1,vHistPca1,'b')
    % stairs(tAxisPca2,vHistPca2,'g')
    xlabel('Month')
    ylabel('Proportion Alive')
    xlim([1 60])
    ylim([0 1])
    title('PCA')

    % subplot(122)
    % stairs(a1,b1,'r')
    % hold on
    % stairs(a2,b2,'b')
    % stairs(a3,b3,'g')
    % xlabel('Month')
    % ylabel('Proportion Alive')
    % xlim([1 60])
    % ylim([0 1])
    % title('PCA (kmplot)')

    [vHistClust0, tAxisClust0] = getVitalHist(vitalStat(groupsClust == 1),month(groupsClust == 1));
    [vHistClust1, tAxisClust1] = getVitalHist(vitalStat(groupsClust ~= 1),month(groupsClust ~= 1));
    % [vHistClust2, tAxisClust2] = getVitalHist(vitalStat(groupsClust == 3),month(groupsClust == 3));

    % figure;
    % [~,~,a1,b1] = kmplot([ month(groupsClust == 1) , vitalStat(groupsClust == 1) ]);
    % [~,~,a2,b2] = kmplot([ month(groupsClust == 2) , vitalStat(groupsClust == 2) ]);
    % [~,~,a3,b3] = kmplot([ month(groupsClust == 3) , vitalStat(groupsClust == 3) ]);

    subplot(122)
    stairs(tAxisClust0,vHistClust0,'r')
    hold on
    stairs(tAxisClust1,vHistClust1,'b')
    % stairs(tAxisClust2,vHistClust2,'g')
    xlabel('Month')
    ylabel('Proportion Alive')
    xlim([1 60])
    ylim([0 1])
    title('Clustering')
    hold off
end
% groupsPca5years = groupsPca;


group1Pca5years = (month <= maxMonth | vitalStat) & groupsPca == 1; %groupsPca == 1;
group2Pca5years = (month <= maxMonth | vitalStat) & groupsPca == 2; %groupsPca == 2;
month1Pca5years = month(group1Pca5years);
month1Pca5years(month1Pca5years > maxMonth & vitalStat(group1Pca5years)) = maxMonth;
month2Pca5years = month(group2Pca5years);
month2Pca5years(month2Pca5years > maxMonth & vitalStat(group2Pca5years)) = maxMonth;

x1 = [ month1Pca5years , vitalStat(group1Pca5years) ];
x2 = [ month2Pca5years , vitalStat(group2Pca5years) ];
figure;
pPca = logrank(x1,x2);

if (plotFlag)
%     figure;
    fSize = 15;
%     title('Kaplan-Meier estimate of survival functions (PCA)','Interpreter','Latex','FontSize',fSize)
    title('')
%     txt1 = strcat('p-Value = ',sprintf('%.4f',pPca));
%     xText = 30; yText = 0.2;
%     text(xText,yText,txt1,'Interpreter','Latex','FontSize',fSize)
    xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
    ylabel('Esimated survival function','Interpreter','Latex','FontSize',fSize)
    le = legend('Group 1','Censord','Group 2')
    set(le,'Interpreter','Latex','FontSize',fSize-2)
    xlim([0 60])
end

group1Clust5years = (month <= maxMonth | vitalStat) & groupsClust == 1; %groupsClust == 1;
group2Clust5years = (month <= maxMonth | vitalStat) & groupsClust ~= 1; %groupsClust == 2;
month1Clust5years = month(group1Clust5years);
month1Clust5years(month1Clust5years > maxMonth & vitalStat(group1Clust5years)) = maxMonth;
month2Clust5years = month(group2Clust5years);
month2Clust5years(month2Clust5years > maxMonth & vitalStat(group2Clust5years)) = maxMonth;

x1 = [ month1Clust5years , vitalStat(group1Clust5years) ];
x2 = [ month2Clust5years , vitalStat(group2Clust5years) ];
figure;
logrank(x1,x2);
if (plotFlag)
%     figure;
%     title('Kaplan-Meier estimate of survival functions (Clustering)','Interpreter','Latex','FontSize',fSize)
%     txt = strcat('p-Value = ',sprintf('%.4f',pClust));
%     txt2 = strcat('p-Value (SPARCoC) = ',sprintf('%.4f',pSPARCoC));
%     xText = 30; yText = 0.2; yText2 = 0.1;
%     text(xText,yText,txt,'Interpreter','Latex','FontSize',fSize)
%     text(xText,yText2,txt2,'Interpreter','Latex','FontSize',fSize)
    title('')
    xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
    ylabel('Esimated survival function','Interpreter','Latex','FontSize',fSize)
    le = legend('Group 1','Censord','Group 2')
    set(le,'Interpreter','Latex','FontSize',fSize-2)
    xlim([0 60])
end
 
% pClustAll(kk) = pClust;
% pPcatAll(kk) = pPca;
votesClust(:,kk) = groupsClust;
VclustAll(:,kk)  = Vclust;

mm
pp
kk

end

% figure;
% plot(knn,pPca,'r')
% hold on
% plot(knn,pClust,'b')
% legend('PCA','Clustering')
% xlabel('K-nn')
% ylabel('P-value')


votesClust(votesClust == 2)       = 0;
votesClust                        = sum(votesClust,2);
% groupsClust                       = ones(N,1);
% groupsClust(votesClust > kIter/2) = 2;
% 
% group1Clust5years = (month <= maxMonth | vitalStat) & groupsClust == 1; %groupsClust == 1;
% group2Clust5years = (month <= maxMonth | vitalStat) & groupsClust ~= 1; %groupsClust == 2;
% month1Clust5years = month(group1Clust5years);
% month1Clust5years(month1Clust5years > maxMonth & vitalStat(group1Clust5years)) = maxMonth;
% month2Clust5years = month(group2Clust5years);
% month2Clust5years(month2Clust5years > maxMonth & vitalStat(group2Clust5years)) = maxMonth;
% 
% x1 = [ month1Clust5years , vitalStat(group1Clust5years) ];
% x2 = [ month2Clust5years , vitalStat(group2Clust5years) ];
% pClust = logrank(x1,x2);


Vclust                  = mean(VclustAll,2);
groupsClust             = ones(N,1);
groupsClust(Vclust > 0) = 2;


group1Clust5years = (month <= maxMonth | vitalStat) & groupsClust == 1; %groupsClust == 1;
group2Clust5years = (month <= maxMonth | vitalStat) & groupsClust ~= 1; %groupsClust == 2;
month1Clust5years = month(group1Clust5years);
month1Clust5years(month1Clust5years > maxMonth & vitalStat(group1Clust5years)) = maxMonth;
month2Clust5years = month(group2Clust5years);
month2Clust5years(month2Clust5years > maxMonth & vitalStat(group2Clust5years)) = maxMonth;

x1 = [ month1Clust5years , vitalStat(group1Clust5years) ];
x2 = [ month2Clust5years , vitalStat(group2Clust5years) ];
pClustAllP(pp) = logrank(x1,x2);
pPcaAllP(pp)   = pPca;

end

pClustAllK(mm)  = mean(pClustAllP);
pPcaAllK(mm)    = mean(pPcaAllP);
pVarClust(mm)   = std(pClustAllP)^2;

end

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
xlim([10 60])
% ylim([0 0.5])


% figure;
% plot(knn,sumD,'k','linewidth',2)
% hold on
% plot(knn(1:end-1),sumD(2:end)+sd_k(2:end))
% hold on
% plot(knn(6),sumD(6),'o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor'...
%     ,0.7*[1 1 1],'MarkerSize',8);
% xlabel('N','Interpreter','Latex','FontSize',fSize)
% ylabel('Gap Statistic, $g(N,8)$','Interpreter','Latex','FontSize',fSize)
% xlim([10 60])
% text(12,3.6,txt,'Interpreter','latex')
% 
% figure;
% plot(knn,pPcaAllK-pClustAllK)
% xlabel('K-nn','Interpreter','Latex','FontSize',fSize)
% ylabel('Difference P-value','Interpreter','Latex','FontSize',fSize)
% title('P-value Difference ($P_{PCA} - P_{Clustering}$)','Interpreter','Latex','FontSize',fSize)
% xlim([10 60])

