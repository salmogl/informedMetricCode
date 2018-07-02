clear all;
close all hidden;
clc; 

addpath('..\');
addpath('kmplot')
addpath('logrank')
addpath('dnafinder-kmplot-cff01a4')
addpath('dnafinder-logrank-5246b53')

dataFolderPathLoad = '';
load(strcat(dataFolderPathLoad,'dataStruct_MSK'));

fSize               = 15;
set(0,'DefaultFigureVisible','off')

data                = dataStruct.genesData;
month               = dataStruct.month;
vitalStat           = dataStruct.vitalStat;

[val,indMaxStd]     = sort(std(data,[],2),'descend');
data                = data(indMaxStd(2:201),:);
[M,N]               = size(data);
rank                = 5;
kClusters           = 21; %10:1:40;

params.rank         = rank;
params.mu           = 1e-1;
params.b            = 1;
params.niter        = 500;
params.stepOpt      = 'constStepLength';
params.debugMode    = true;
params.threshDiff   = 5e-5;

nIter               = 1;
pPcaAllK            = zeros(1,length(kClusters));
pInfAllK            = zeros(1,length(kClusters));
pVarInfAllK         = zeros(1,length(kClusters));

for kk = 1:length(kClusters)

    pInfAllp            = zeros(1,nIter);
    pPcaAllp            = zeros(1,nIter);

    for ii = 1 : nIter

        class               = kmeans(data,kClusters(kk),'MaxIter',1e6,'OnlinePhase','on','replicates',1e5);
        [~,indSort]         = sort(class);

        [~,S0,Upca]         = svd(cov(data'));
        Upca                = Upca(:,1:rank);
        S                   = S0(1:rank,1:rank);
        covPca              = Upca*S*Upca';
        covSam              = cov(data');
        Uinf                = getPdSubProj(data,class,params);
        covInf              = Uinf*S*Uinf';

        %%
        dataCentered        = bsxfun(@minus,data,mean(data,2));
        Vpca                = (Upca*S^(-1/2))'*dataCentered;
        Vinf                = (Uinf*S^(-1/2))'*dataCentered;

        %%
        groupsPca           = kmeans(Vpca',2,'MaxIter',1e6,'OnlinePhase','on','replicates',500);
        [x1_pca,x2_pca]     = getRiskGroups(month,vitalStat,groupsPca);
        pPcaAllp(ii)        = logrankWrap(x2_pca,x1_pca);

        groupsInf           = kmeans(Vinf',2,'MaxIter',1e6,'OnlinePhase','on','replicates',500);
        [x1_inf,x2_inf]     = getRiskGroups(month,vitalStat,groupsInf);
        pInfAllp(ii)        = logrankWrap(x2_inf,x1_inf);

    end

    pPcaAllK(kk)         = mean(pPcaAllp);
    pInfAllK(kk)         = mean(pInfAllp);
    pVarInfAllK(kk)      = var(pInfAllp);

end

%% Figures
set(0,'DefaultFigureVisible','on')

% Data
% genes are ordered according to their variance
HMObject = HeatMap(data);
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)

% genes are ordered according to clusters obtained by k-means
HMObject = HeatMap(data(indSort,:));
addYLabel(HMObject, 'Genes','FontSize', fSize)
addXLabel(HMObject, 'Subjects','FontSize', fSize)


% Eigenvalues decay
dataSet             = 'CANDF';
loadPath            = addPathLoad(dataSet);
load(loadPath)
dataCANDF           = dataStruct.genesData;
[val,indMaxStd]     = sort(std(dataCANDF,[],2),'descend');
dataCANDF           = dataCANDF(indMaxStd(2:201),:);
[~,S0CANDF,~]       = svd(cov(dataCANDF'));
eigsMSK             = sum(S0)/max(S0(:));
eigsCAN             = sum(S0CANDF)/max(S0CANDF(:));

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

% Kaplan-Meier survival plots
figure;
pPca                = logrankWrap(x2_pca,x1_pca);
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Estimated survival function','Interpreter','Latex','FontSize',fSize)
le                  = legend('Group 1','Censord','Group 2');
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])

figure;
pInf                = logrankWrap(x2_inf,x1_inf);
title('')
xlabel('Time [Month]','Interpreter','Latex','FontSize',fSize)
ylabel('Estimated survival function','Interpreter','Latex','FontSize',fSize)
le                  = legend('Group 1','Censord','Group 2');
set(le,'Interpreter','Latex','FontSize',fSize-2)
xlim([0 60])


% P-Value Vs. number of clusters
figure;
plot(kClusters,log(pPcaAllK));
hold on
plot(kClusters,log(pInfAllK));
ylabel('p-Value [dB]','fontsize',fSize)
xlabel('K [Clusters]','fontsize',fSize)
xlim([10 30])
le = legend('$\hat{p}$','$\tilde{p}$');
set(le,'Interpreter','Latex','FontSize',fSize)

% Comparison of covariance matrices
figure;
imagesc(covSam(indSort,indSort))
% title('Sample Covariance')

figure;
imagesc(covPca(indSort,indSort));
% title('PCA')

figure;
imagesc(covInf(indSort,indSort))
% title('Informed PCA')