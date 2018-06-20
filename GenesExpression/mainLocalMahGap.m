clear all;
close all hidden;
clc;

dataSet  = 'MSK';
loadPath = addPathLoad(dataSet);
load(loadPath);
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

color = jet(max(month)+1);
colorsVital = color(month,:);
colorsVital(vitalStat == 1,:) = repmat(color(max(month)+1,:),length(find(vitalStat)), 1);
colorsStage = dataStruct.stage; 

[val,deadOrder] = sort(month(~vitalStat));
colorsVital0 = jet(length(month(~vitalStat)));
[~,kLargest] = sort(deadOrder);
colorsVital0 = colorsVital0(kLargest,:);

deadIdx = 1:N;
deadIdx = deadIdx(~vitalStat);
 

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

knn         = [10 14:2:60]; %[45:5:75 82];%[10 15 20 25 30 35 40 45]; % 20 is the optimal
pClustAllP  = zeros(1,pIter);
pPcaAllP    = zeros(1,pIter);
pClustAllK  = zeros(1,length(knn));
pVarClust   = zeros(1,length(knn));
pPcaAllK    = zeros(1,length(knn));
sumD        = zeros(1,length(knn));
sd_k        = zeros(1,length(knn));

argsPca.knn            = knn;
argsPca.replicates     = 10; %500
argsPca.nClust         = nClust;
argsPca.rank           = rank;
argsPca.informedFlag   = false;
argsInf                = argsPca;
argsInf.informedFlag   = true;

argsDiff.eps        = 1;
argsDiff.eigsnum    = 80;

parfor mm = 1:length(knn)
    
    knn(mm)
    argsPcaCur          = argsPca;
    argsPcaCur.knn      = knn(mm);
    argsInfCur          = argsInf;
    argsInfCur.knn      = knn(mm);
    mDistPca            = PCA_ReconstructionGap(data,argsPcaCur);
  
    embeddingPca        = diffusionMap(mDistPca,argsDiff);
    [~,~,Vpca]          = svd(embeddingPca');
    groupsPca           = ones(N,1);
    Vpca                = Vpca(:,1)*Vpca(1,1);
    groupsPca(Vpca > 0) = 2;
    
    [x1_pca,x2_pca]     = getRiskGroups(month,vitalStat,groupsPca);
    pPca                = logrank(x1_pca,x2_pca);

    if (plotFlag)
        figure;
        logrank(x1_pca,x2_pca);
        fSize = 15;
        title('')
        xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
        ylabel('Esimated survival function','Interpreter','Latex','FontSize',fSize)
        le = legend('Group 1','Censord','Group 2')
        set(le,'Interpreter','Latex','FontSize',fSize-2)
        xlim([0 60])
    end
    
    
    
for pp = 1:pIter
    
    votesClust          = zeros(N,kIter);
    VclustAll           = zeros(N,kIter);
    
for kk = 1:kIter
    
[mDistClust,sumD(mm)] = PCA_ReconstructionGap(data,argsInfCur);
% [~,mDistClust]        = PCA_Reconstruction(data,args);
embeddingClust        = diffusionMap(mDistClust,argsDiff);
[~,~,Vclust]          = svd(embeddingClust');
% votesClust(:,kk)      = groupsClust;
Vclust                = Vclust(:,1)*sign(Vclust(1,1));
VclustAll(:,kk)       = Vclust;


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


mm
pp
kk

end


votesClust(votesClust == 2)   = 0;
votesClust                    = sum(votesClust,2);
Vclust                        = mean(VclustAll,2);
groupsClust                   = ones(N,1);
groupsClust(Vclust > 0)       = 2;
[x1_clust,x2_clust]           = getRiskGroups(month,vitalStat,groupsClust);
%pClustAllP(pp)                = logrank(x1_clust,x2_clust);

end

%pClustAllK(mm)  = mean(PClustAllP);
pPcaAllK(mm)    = pPca;
%pVarClust(mm)   = std(pClustAllP)^2;

end

fSize = 18;
std0 = sqrt(pVarClust);
figure;
% plot(x,y,'r.-');
h(1) = plot(knn,pPcaAllK,'r.-')
hold on
h(2) = patch([knn fliplr(knn)],[pClustAllK+std0 fliplr(pClustAllK-std0)],0.8*[1 1 1]);
h(3) = plot(knn,pClustAllK,'b.-')

le = legend(h([1 3]),'$p_{LM}$','$p_{ILM}$')

set(le,'Interpreter','Latex','FontSize',fSize-3)
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('P-value','Interpreter','Latex','FontSize',fSize)
xlim([10 60])
ylim([0 0.5])

%%
figure;
plot(knn,sumD,'k','linewidth',2)
hold on
% plot(knn(1:end-1),sumD(2:end)+sd_k(2:end))
% hold on
plot(knn(12),sumD(12),'o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor'...
    ,0.7*[1 1 1],'MarkerSize',8);
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('Gap Statistic, $g(N,8)$','Interpreter','Latex','FontSize',fSize)
xlim([10 60])
% text(12,3.6,txt,'Interpreter','latex')
% 
% figure;
% plot(knn,pPcaAllK-pClustAllK)
% xlabel('K-nn','Interpreter','Latex','FontSize',fSize)
% ylabel('Difference P-value','Interpreter','Latex','FontSize',fSize)
% title('P-value Difference ($P_{PCA} - P_{Clustering}$)','Interpreter','Latex','FontSize',fSize)
% xlim([10 60])

