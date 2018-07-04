% Source code for the paper "Mahalanonbis Distance Informed by Clustering" 
% by Almog Lahav, Ronen Talmon and Yuval Kluger.
%==========================================================================
% Section 7 (in the paper) - "Application to Gene Expression Data" 
% Figures  : 15a, 15b, 16, 17
% Comments : 
% To obtain Fig. 16 and Fig. 17 - set:
%                   knn 	= 10 : 5 :60; 
%                   pIter 	= 20;
% Note that runing with this setting will take a long time.  

clear all;
close all hidden;
clc; 

addpath('..\')
addpath('dnafinder-kmplot-cff01a4')
addpath('dnafinder-logrank-5246b53')

dataFolderPathLoad = 'C:\Users\myPath\';
load(strcat(dataFolderPathLoad,'dataStruct_CANDF'));

fSize       = 15;
plotFlag    = false;

if (plotFlag)
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

pSPARCoC        = 0.0096; % p-Value obtained by SPARCoC method
data            = dataStruct.genesData;
month           = dataStruct.month;
vitalStat       = dataStruct.vitalStat;
[~,indMaxStd]   = sort(std(data,[],2),'descend');
data            = data(indMaxStd(2:201),:);
[M,N]           = size(data);
rank            = 6;
nClust          = 7;

[class,centroids] = kmeans(data,nClust,'MaxIter',1e6,'OnlinePhase','on');
[~,indSort]       = sort(class);

HMObject = HeatMap(data(indSort,:));
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)

HMObject = HeatMap(data);
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)


%%
kIter           = 1;                    % PC of people is avaraged over kIter iterations
pIter           = 1;                    % p-value is avaraged over pIter iterations. Figures 16-17 in the paper: pIter = 20;
knn             = 34;                   % neighborhood size. Figures 16-17 in the paper: knn  = 10 : 5 :60;
pClustAllP      = zeros(1,pIter);
pPcaAllP        = zeros(1,pIter);
pClustAllK      = zeros(1,length(knn));
pVarClust       = zeros(1,length(knn));
pPcaAllK        = zeros(1,length(knn));

args.knn        = knn;
args.replicates = 300;
args.nClust     = nClust;
args.rank       = rank;
maxMonth        = 60;

for mm = 1:length(knn)
    
    args.knn    = knn(mm);

    for pp = 1:pIter

        VclustAll   = zeros(N,kIter);

        for kk = 1:kIter
            
            [mDistPca,mDistClust]   = localMahalanobis(data,args);
            eigsnum                 = 80;
            eps                     = 1;
            sigma                   = eps * median(mDistClust(:));
            affClust                = exp(-mDistClust.^2/(2*sigma^2)); 
            affClust                = stochastic(affClust);
            [eigvecs, eigvals]      = eigs(affClust, eigsnum);
            embeddingClust          = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);
            [~,~,Vclust]            = svd(embeddingClust');
            VclustAll(:,kk)         = Vclust(:,1)*sign(Vclust(1,1));


        end

        sigma                   = eps * median(mDistPca(:));
        affPca                  = exp(-mDistPca.^2/(2*sigma^2)); 
        affPca                  = stochastic(affPca);
        [eigvecs, eigvals]      = eigs(affPca, eigsnum);
        embeddingPca            = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

        deadIdx                 = 1:N;
        deadIdx                 = deadIdx(~vitalStat);
        [~,~,Vpca]              = svd(embeddingPca');
        groupsPca               = ones(N,1);
        Vpca                    = Vpca(:,1)*Vpca(1,1);
        groupsPca(Vpca > 0)     = 2;
        group1Pca5years         = (month <= maxMonth | vitalStat) & groupsPca == 1;
        group2Pca5years         = (month <= maxMonth | vitalStat) & groupsPca == 2;
        month1Pca5years         = month(group1Pca5years);
        month1Pca5years(month1Pca5years > maxMonth & vitalStat(group1Pca5years)) = maxMonth;
        month2Pca5years         = month(group2Pca5years);
        month2Pca5years(month2Pca5years > maxMonth & vitalStat(group2Pca5years)) = maxMonth;
        x1Lm                    = [ month1Pca5years , vitalStat(group1Pca5years) ];
        x2Lm                    = [ month2Pca5years , vitalStat(group2Pca5years) ];
        pPcaAllP(pp)            = logrankWrap(x1Lm,x2Lm);

        Vclust                  = mean(VclustAll,2);
        groupsClust             = ones(N,1);
        groupsClust(Vclust > 0) = 2;
        group1Clust5years       = (month <= maxMonth | vitalStat) & groupsClust == 1; %groupsClust == 1;
        group2Clust5years       = (month <= maxMonth | vitalStat) & groupsClust ~= 1; %groupsClust == 2;
        month1Clust5years       = month(group1Clust5years);
        month1Clust5years(month1Clust5years > maxMonth & vitalStat(group1Clust5years)) = maxMonth;
        month2Clust5years       = month(group2Clust5years);
        month2Clust5years(month2Clust5years > maxMonth & vitalStat(group2Clust5years)) = maxMonth;
        x1Ilm                   = [ month1Clust5years , vitalStat(group1Clust5years) ];
        x2Ilm                   = [ month2Clust5years , vitalStat(group2Clust5years) ];
        pClustAllP(pp)          = logrankWrap(x1Ilm,x2Ilm);

    end

    pClustAllK(mm)  = mean(pClustAllP);
    pPcaAllK(mm)    = mean(pPcaAllP);
    pVarClust(mm)   = std(pClustAllP)^2;

end

set(0,'DefaultFigureVisible','on')

figure;
logrankWrap(x1Lm,x2Lm);
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Esimated survival function','Interpreter','Latex','FontSize',fSize)
le = legend('Group 1','Censord','Group 2');
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])

figure;
logrankWrap(x1Ilm,x2Ilm);
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Esimated survival function','Interpreter','Latex','FontSize',fSize)
le = legend('Group 1','Censord','Group 2');
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])

%%
fSize   = 18;
std0    = sqrt(pVarClust);

figure;
h(1)    = plot(knn,pPcaAllK,'r.-');
hold on
h(2)    = patch([knn fliplr(knn)],[pClustAllK+std0 fliplr(pClustAllK-std0)],0.8*[1 1 1]);
h(3)    = plot(knn,pClustAllK,'b.-');
le      = legend(h([1 3]),'$p_{LM}$','$p_{ILM}$');
set(le,'Interpreter','Latex','FontSize',fSize-3)
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('P-value','Interpreter','Latex','FontSize',fSize)
xlim([10 60])

%%
argsGap.nClust     = nClust;
argsGap.replicates = 100;
sumD               = zeros(1,length(knn));

for ii = 1:length(knn)
    argsGap.knn     = knn(ii);
    sumD(ii)        = gapStat(data,argsGap);
end

figure;
plot(knn,sumD,'k','linewidth',2)
xlabel('N','Interpreter','Latex','FontSize',fSize)
ylabel('Gap Statistic, $g(N,7)$','Interpreter','Latex','FontSize',fSize)
xlim([10 60])