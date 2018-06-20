clear all;
close all hidden;
clc; 
dataSet             = 'MSK';
loadPath            = addPathLoad(dataSet);
load(loadPath);

plotFlag            = true;
if (plotFlag)
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

data                = dataStruct.genesData;
month               = dataStruct.month;
vitalStat           = dataStruct.vitalStat;

[val,indMaxStd]     = sort(std(data,[],2),'descend');
data                = data(indMaxStd(2:201),:);
[M,N]               = size(data);
rank                = 5;
kClusters           = 21;%10:1:30;

for kk = 1:length(kClusters)

class               = kmeans(data,kClusters(kk),'MaxIter',1e6,'OnlinePhase','on','replicates',1e2);
[~,indSort]         = sort(class);

figure;
fSize = 15;
HMObject = HeatMap(data(indSort,:));
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)

figure;
HMObject = HeatMap(data);
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)



[~,S0,Upca]         = svd(cov(data'));
Upca                = Upca(:,1:rank);
S                   = S0(1:rank,1:rank);
covPca              = Upca*S*Upca';
covSam              = cov(data');


params.rank         = rank;
params.mu           = 1e-1;
params.b            = 1;
params.niter        = 500;
params.stepOpt      = 'constStepLength';
params.debugMode    = true;
params.threshDiff   = 5e-5;
Uinf                = getPdSubProj(data,class,params);
covInf              = Uinf*S*Uinf';

%%
figure;
imagesc(covPca(indSort,indSort));
% title('PCA')

figure;
imagesc(covSam(indSort,indSort))
% title('Sample Covariance')

figure;
imagesc(covInf(indSort,indSort))
% title('Informed PCA')

%%
dataCentered = bsxfun(@minus,data,mean(data,2));
Vpca         = (Upca*S^(-1/2))'*dataCentered;
Vinf         = (Uinf*S^(-1/2))'*dataCentered;
figure;
subplot(1,2,1)
plot3dData(Vpca(1:3,:),vitalStat,'Embedding Colored By Vital Status (Dead) (PCA)')

subplot(1,2,2)
plot3dData(Vinf(1:3,:),vitalStat,'Embedding Colored By Vital Status (Dead) (Informed)')

%%
% groupsPca           = ones(N,1);
% Vpca1               = Vpca(1,:)*Vpca(1,1);
% groupsPca(Vpca1 > 0)= 2;
groupsPca           = kmeans(Vpca',2,'MaxIter',1e6,'OnlinePhase','on','replicates',500);
%%
[x1_pca,x2_pca]     = getRiskGroups(month,vitalStat,groupsPca);
figure;
pPca                = logrank(x2_pca,x1_pca);
fSize               = 15;
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Estimated survival function','Interpreter','Latex','FontSize',fSize)
le                  = legend('Group 1','Censord','Group 2')
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])

%%
% groupsInf           = ones(N,1);
% Vinf1               = Vinf(1,:)*Vinf(1,1);
% groupsInf(Vinf1 > 0)= 2;tsneVc

groupsInf           = kmeans(Vinf',2,'MaxIter',1e6,'OnlinePhase','on','replicates',500);
%%
[x1_inf,x2_inf]     = getRiskGroups(month,vitalStat,groupsInf);
figure;
pInf                = logrank(x2_inf,x1_inf);
fSize               = 15;
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Estimated survival function','Interpreter','Latex','FontSize',fSize)
le                  = legend('Group 1','Censord','Group 2')
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])

pPcaAllK(kk) = pPca;
pInfAllK(kk) = pInf;

end
%%
% no_dims         = 3;
% initial_dims    = rank;
% perplexity      = 20;
% clustLabels     = vitalStat;
% figure;
% tsneVc          = tsne_inf(data', clustLabels, no_dims, initial_dims, perplexity);
% tsneRaw         = tsne(Vinf', clustLabels, no_dims, initial_dims, perplexity);
%%
dataSet             = 'CANDF';
loadPath            = addPathLoad(dataSet);
load(loadPath)
dataCANDF           = dataStruct.genesData;
[val,indMaxStd]     = sort(std(dataCANDF,[],2),'descend');
dataCANDF           = dataCANDF(indMaxStd(2:201),:);
[~,S0CANDF,~]       = svd(cov(dataCANDF'));
%%
eigsMSK         = sum(S0)/max(S0(:));
eigsCAN         = sum(S0CANDF)/max(S0CANDF(:));

figure;
plot(eigsMSK)
hold on
plot(eigsCAN)
le = legend('MSK','CANDF');
set(le,'Interpreter','Latex','FontSize',fSize)
plot(ones(1,length(eigsCAN))*0.15,'k--')
plot(ones(1,length(eigsCAN))*6,linspace(0,eigsCAN(6),length(eigsCAN)),'r--')
plot(ones(1,length(eigsMSK))*5,linspace(0,eigsMSK(5),length(eigsMSK)),'b--')
xlim([1 20])
ylabel('Normalized eigenvalues ($\lambda_i$)','Interpreter','Latex','FontSize',fSize)
xlabel('$i$','Interpreter','Latex','FontSize',fSize)
%%
fSize = 15;
load('../MatFiles/pValVsK_10_1_40_var.mat')
% pPcaAllK_70_5_100   = pPcaAllK;
% pInfAllK_70_5_100   = pInfAllK;
% kClusters_70_5_100  = kClusters;
% load('../MatFiles/pValVsK_1.mat')
% pPcaAllK            = [pPcaAllK pPcaAllK_70_5_100];
% pInfAllK            = [pInfAllK pInfAllK_70_5_100];
% kClusters           = [kClusters kClusters_70_5_100];
std0  = pVarInfAllK;
kClusters = kClusters(1:length(pPcaAllK))
figure;
% plot(x,y,'r.-');
h(1) = plot(kClusters,pPcaAllK,'r.-')
hold on
h(2) = patch([kClusters fliplr(kClusters)],[pInfAllK+std0 fliplr(pInfAllK-std0)],0.8*[1 1 1]);
h(3) = plot(kClusters,pInfAllK,'b.-')

le = legend(h([1 3]),'$\hat{p}$','$\tilde{p}$')

set(le,'Interpreter','Latex','FontSize',fSize-3)

figure;
plot(kClusters,log(pPcaAllK));
hold on
plot(kClusters,log(pInfAllK));
ylabel('p-Value [dB]','fontsize',15)
xlabel('K [Clusters]','fontsize',15)
xlim([10 30])
le = legend('$\hat{p}$','$\tilde{p}$');
set(le,'Interpreter','Latex','FontSize',15)