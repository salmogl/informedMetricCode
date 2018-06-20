clc;
clear all;
close all hidden;
close all;
dbstop if error

addpath('C:\Users\salmogl\Google Drive\Master\MATLAB\CoupeledGeomrtry\3DQuest\3DQuest\Questionnaire')

% addpath('../3DQuest/3DQuest/Questionnaire')

set(0,'defaulttextinterpreter','latex')

nTrial      = 1;
nColsTest   = 100;%[20:20:200 250:50:700 800:100:2000];%[50 75 100 200 400 800];
plotFlag    = false;

if (plotFlag)
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

VkClassRows = 18;

for iter = 1:length(VkClassRows);

errAllNorms = zeros(4,3);
VkClassRows(iter)

% Data Generation
args.option = 5;
args.nRows = 900;
args.nCols = nColsTest;
args.nColsTest = nColsTest;
args.kClassRows = 18;
args.sigma = 0.5;
args.sigmaNoise = 8;
args.randFlag = false;
args.sigmaNoise = 0;

% Affifnity parameters
paramsAff.eps        = 1;
paramsAff.knn        = nColsTest;
paramsAff.metric     = 'euclidean'; % 'cosine_similarity'/'euclidean'
paramsAff.thresh     = 0.05;

% Projected gradient parameters
paramsGrad.rank      = args.kClassRows;
paramsGrad.mu        = 3e-4;
paramsGrad.niter     = 200;
paramsGrad.stepOpt   = 'diminish';
paramsGrad.b         = 1;
paramsGrad.debugMode = false;
paramsGrad.mu        = 1e-1;
paramsGrad.stepOpt   = 'constStepLength';
paramsGrad.threshDiff = 0;

kClassRows = args.kClassRows;

for trial = 1:nTrial
    
% Simulation Parameters
embedded   = false;
coupling = false;
dataTypa = 'noPerm';
plotFlagKmeans = plotFlag;
useFit = false;

[dataTest,bordersRows,data,origData,covariance,columnsClass,columnsClassTest] = getData(args);

data = data + args.sigmaNoise*rand(size(data));
dataTest = dataTest + args.sigmaNoise*rand(size(dataTest));

%%%%%% Algorithm Start %%%%%%
data = data - repmat(mean(data,2),1,size(data,2));
dataTest = dataTest - repmat(mean(dataTest,2),1,size(dataTest,2));

[realPc,sReal,v] = svd(covariance);
realPcNorm = realPc/sqrt(sReal);
[sizeRows,sizeCols] = size(dataTest);

switch(dataTypa)
    case 'noPerm'
        shuffleRows = 1:sizeRows;
        shuffleCols = 1:sizeCols; 
    case 'randPerm'
        shuffleRows = randperm(sizeRows);
        shuffleCols = randperm(sizeCols);
    case 'constPerm'
%         load('Animation/permutations1');
%         shuffleRows = randperm(sizeRows);
%         shuffleCols = randperm(sizeCols); 
end

shuffleCols = columnsClassTest;

% data = dataOrig(shuffleRows,shuffleCols);

%%% Embedding computation using euclidean/cosine distance 
aff = CalcInitAff(dataTest,paramsAff);

aff = stochastic(aff);
eigsnum = 4;
[eigvecs, eigvals] = eigs(aff, eigsnum);
embedding = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

%%% Embedding computation using Mahalanobis distance 
% kClassRows = 7;
%[class,V,centroids,S] = KmeansSvd(dataTest',kClassRows,useFit,plotFlagKmeans,'');
tic
[class,~,D] = kmeans(dataTest,VkClassRows(iter),'OnlinePhase','on','Replicates',500);
toc
[V,S,U] = svd(cov(dataTest'));

% class = kmeans(dataTest,kClassRows);
dataTestProjectedMahalanobis = (V(:,1:kClassRows-1)/sqrt(S(1:kClassRows-1,1:kClassRows-1)))'*dataTest;

aff = CalcInitAff(dataTestProjectedMahalanobis,paramsAff);
aff = stochastic(aff);
eigsnum = 4;
[eigvecs, eigvals] = eigs(aff, eigsnum);
embeddingMahalanobis = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

%%% Embedding computation using Mahalanobis-Like distance (clusters using k-means)
classNum = length(unique(class));
Classifier = zeros(size(class,1), classNum);
for i=1:classNum
    Classifier(class == i,i) = 1;
end

% normalization according to group size
ClassifierNorm = Classifier./(repmat(sqrt(sum(Classifier)),[size(class,1) 1]));
Q = [ones(sizeRows,1)/sqrt(sizeRows) V(:,1:classNum-1)];
T = ((ClassifierNorm'*ClassifierNorm)^-1)*ClassifierNorm'*Q;
T = T(:,2:end)/sqrt(S(1:classNum-1,1:classNum-1));
dataTestProjectedClass = (ClassifierNorm*T)'*dataTest;

dataTestProjectedClass0 = (ClassifierNorm)'*dataTest;
[V0,S0,U0] = svd(cov(dataTestProjectedClass0'));


aff = CalcInitAff(dataTestProjectedClass,paramsAff);
aff = stochastic(aff);
eigsnum = 4;
[eigvecs, eigvals] = eigs(aff, eigsnum);
embeddingKmeans = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

%%% Embedding computation using Mahalanobis-Like distance (known clusters)

classNumReal    = kClassRows;
classSize       = floor(sizeRows/classNumReal);
Classifier      = zeros(sizeRows, classNumReal);
[vals,index]    = sort(shuffleRows);
classReal       = zeros(sizeRows,1);
classRealRows   = zeros(sizeRows,1);
for i=1:classNumReal
        Classifier(index(bordersRows(i)+1:bordersRows(i+1)),i) = 1;
        classRealRows(index(bordersRows(i)+1:bordersRows(i+1))) = i;
end

% normalization according to group size
ClassifierNormReal = Classifier./(repmat(sqrt(sum(Classifier)),[sizeRows 1]));
Q = [ones(sizeRows,1)/sqrt(sizeRows) V(:,1:classNumReal-1)];
TReal = ((ClassifierNormReal'*ClassifierNormReal)^-1)*ClassifierNormReal'*Q;
TReal = TReal(:,2:end)/sqrt(S(1:classNumReal-1,1:classNumReal-1));
dataTestProjectedRealClass = (ClassifierNormReal*TReal)'*dataTest;
% dataTestProjectedRealClass = (ClassifierNormReal)'*dataTest;
dataTestProjectedRealClass1 = (ClassifierNormReal)'*dataTest;
[V1,S1,U1] = svd(cov(dataTestProjectedRealClass1'));

aff = CalcInitAff(dataTestProjectedRealClass,paramsAff);
aff = stochastic(aff);
eigsnum = 4;
[eigvecs, eigvals] = eigs(aff, eigsnum);
embeddingRealClusters = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

V_real = realPc(:,1:kClassRows);
V_pca  = V(:,1:kClassRows);
% V_inf  = ClassifierNorm*ClassifierNorm'*V_pca;
% V_infR = ClassifierNormReal*ClassifierNormReal'*V_pca;
%%

V_inf    = getPdSubProj(dataTest,class,paramsGrad);
% V_inf    = V_inf./repmat(sqrt(sum(V_inf.^2,1)),sizeRows,1);

V_infR   = getPdSubProj(dataTest,classRealRows,paramsGrad);
% V_infR   = V_infR./repmat(sqrt(sum(V_infR.^2,1)),sizeRows,1);

S_Pca    = S(1:kClassRows,1:kClassRows);
pc       = V(:,1:kClassRows-1)/sqrt(S(1:kClassRows-1,1:kClassRows-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plotFlag)

figure;
imagesc(dataTest), colormap jet, axis on

figure;
pointSize = 15;
scatter3(embedding(:,1),embedding(:,2),embedding(:,3),pointSize,shuffleCols,'filed')
xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
colormap winter

figure;
scatter3(embeddingMahalanobis(:,1),embeddingMahalanobis(:,2),embeddingMahalanobis(:,3),pointSize,shuffleCols,'filed')
xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
colormap winter

figure;
scatter3(embeddingKmeans(:,1),embeddingKmeans(:,2),embeddingKmeans(:,3),pointSize,shuffleCols,'filed')
xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
colormap winter

figure;
scatter3(embeddingRealClusters(:,1),embeddingRealClusters(:,2),embeddingRealClusters(:,3),pointSize,shuffleCols,'filed')
xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
colormap winter


for i = 1:1
    figure;
    temp = V_pca(:,i)*V_pca(:,i)';
    plot(temp(1,:))
    hold on
    temp = V_inf(:,i)*V_inf(:,i)';
    plot(temp(1,:))    
    temp = V_infR(:,i)*V_infR(:,i)';
    plot(temp(1,:))    
    temp = V_real(:,i)*V_real(:,i)';
    plot(temp(1,:))
     title(strcat('PC #',num2str(i)))
    legend('PCA','PCA With Clustering','PCA With Known Clustering','Real PC')

end

covPCA   = V_pca*S_Pca*V_pca';
% S_inf    = diag(std(dataTest'*V_inf).^2);
covClass = V_inf*S_Pca*V_inf';

%%
plotCovariance(covPCA,'$\widehat{\Sigma}_{ij}$');
plotCovariance(covClass,'$\widetilde{\Sigma}_{ij}$');
plotCovariance(covariance,'$\Sigma_{ij}$');

%%
jrow = 900;
figure;
plot(1:sizeRows,covClass(:,jrow),'linewidth',2)
hold on;
plot(1:sizeRows,covariance(:,jrow),'linewidth',2)
plot(1:sizeRows,covPCA(:,jrow),'k')

xlabel('$i$','fontsize',20)
h=legend('$\widetilde{\Sigma}_{i900}$','$\Sigma_{i900}$','$\widehat{\Sigma}_{i900}$')
set(h,'Interpreter','latex')
set(h,'Position',[0.2 0.7 0.2 0.2])
set(h,'fontsize',15)
ylim([-1 29])


dataProjectedClass              = (ClassifierNorm*T)'*dataTest;
dataProjectedClassNoT           = (ClassifierNorm)'*dataTest;
dataProjectedRealClass          = (ClassifierNormReal*TReal)'*dataTest;
dataProjectedRealClassNoT       = (ClassifierNormReal)'*dataTest;
dataProjectedMahalanobis        = (V(:,1:kClassRows-1)/sqrt(S(1:kClassRows-1,1:kClassRows-1)))'*dataTest;
dataProjectedMahalanobisNoN     = V(:,1:kClassRows-1)'*dataTest;
dataProjectedRealMahalanobis    = (realPcNorm(:,1:kClassRows-1))'*dataTest;


classDist                       = pdist(dataProjectedClass');
classDistNoT                    = pdist(dataProjectedClassNoT');
classRealDist                   = pdist(dataProjectedRealClass');
classRealDistNoT                = pdist(dataProjectedRealClassNoT');
mahalanobisDist                 = pdist(dataProjectedMahalanobis');
mahalanobisDistNoN              = pdist(dataProjectedMahalanobisNoN');
mahalanobisRealDist             = pdist(dataProjectedRealMahalanobis');


M = corrcoef([mahalanobisRealDist' classDist']);
classCor = M(2,1);
M = corrcoef([mahalanobisRealDist' classDistNoT']);
classCorNoT = M(2,1);
M = corrcoef([mahalanobisRealDist' classRealDist']);
classRealCor = M(2,1);
M = corrcoef([mahalanobisRealDist' classRealDistNoT']);
classRealCorNoT = M(2,1);
M = corrcoef([mahalanobisRealDist' mahalanobisDist']);
mahalanobisCor = M(2,1);
M = corrcoef([mahalanobisRealDist' mahalanobisDistNoN']);
mahalanobisCorNoN = M(2,1);

classDist = classDist/max(classDist);
classDistNoT = classDistNoT/max(classDistNoT);
classRealDist = classRealDist/max(classRealDist);
classRealDistNoT = classRealDistNoT/max(classRealDistNoT);
mahalanobisDist = mahalanobisDist/max(mahalanobisDist);
mahalanobisDistNoN = mahalanobisDistNoN/max(mahalanobisDistNoN);
mahalanobisRealDist = mahalanobisRealDist/max(mahalanobisRealDist);


pointSize = 10;
figure;

subplot(3,1,1)
scatter(mahalanobisRealDist,mahalanobisDist,pointSize,'filled')
h=title('Mahalanobis Distance: Sample Covariance Vs. Real Covariance');
set(h,'Interpreter','latex','FontSize',15)
xlabel('Mahalanobis With Real Covariance','Interpreter','latex','FontSize',12)
% ylabel('$\|f\left(c_{i}\right)-f\left(c_{j}\right)\|^{2}$','Interpreter','latex','FontSize',15)
ylabel('Mahalanobis Dist.','Interpreter','latex','FontSize',12)
t=text(0.2,0.8,strcat('Corr = ',sprintf('%.2f',mahalanobisCor)));
t.FontSize = 12;

subplot(3,1,2)
scatter(mahalanobisRealDist,classDist,pointSize,'filled')
h=title('Clustering Distance Vs. Real Mahalanobis');
set(h,'Interpreter','latex','FontSize',15)
xlabel('Mahalanobis With Real Covariance','Interpreter','latex','FontSize',12)
% ylabel('$\|f\left(c_{i}\right)-f\left(c_{j}\right)\|^{2}$','Interpreter','latex','FontSize',15)
ylabel('Clustering Dist.','Interpreter','latex','FontSize',12)
t=text(0.2,0.8,strcat('Corr = ',sprintf('%.2f',classCor)));
t.FontSize = 12;


subplot(3,1,3)
scatter(mahalanobisRealDist,classRealDist,pointSize,'filled')
h=title('Known Clustering Distance Vs. Real Mahalanobis');
set(h,'Interpreter','latex','FontSize',15)
xlabel('Mahalanobis With Real Covariance','Interpreter','latex','FontSize',12)
ylabel('Known Clust. Dist.','Interpreter','latex','FontSize',12)
t=text(0.2,0.8,strcat('Corr = ',sprintf('%.2f',classRealCor)));
t.FontSize = 12;

end

errNorm2Class           = norm(V_inf*V_inf'-V_real*V_real');
errNorm2RealClass       = norm(V_infR*V_infR'-V_real*V_real');
errNorm2Mahalanobis     = norm(V_pca*V_pca'-V_real*V_real');


errNorm1Class           = norm(V_inf*V_inf'-V_real*V_real',1);
errNorm1RealClass       = norm(V_infR*V_infR'-V_real*V_real',1);
errNorm1Mahalanobis     = norm(V_pca*V_pca'-V_real*V_real',1);

errNormInfClass         = norm(V_inf*V_inf'-V_real*V_real',inf);
errNormInfRealClass     = norm(V_infR*V_infR'-V_real*V_real',inf);
errNormInfMahalanobis   = norm(V_pca*V_pca'-V_real*V_real',inf);

errNormFroClass         = norm(V_inf*V_inf'-V_real*V_real','fro');
errNormFroRealClass     = norm(V_infR*V_infR'-V_real*V_real','fro');
errNormFroMahalanobis   = norm(V_pca*V_pca'-V_real*V_real','fro');


errAllNorms = errAllNorms +[errNorm1Class, errNorm1RealClass,errNorm1Mahalanobis;...
    errNorm2Class, errNorm2RealClass, errNorm2Mahalanobis;...
    errNormFroClass, errNormFroRealClass, errNormFroMahalanobis;...
    errNormInfClass, errNormInfRealClass, errNormInfMahalanobis;
    ];


end

errAllNormsMean{iter,1} = errAllNorms/nTrial;

errAllNormsMean{iter,2} = strcat('Number of clusters: ',num2str(VkClassRows(iter)));

end
set(0,'DefaultFigureVisible','on')

%%
errArray = zeros(length(VkClassRows),3,4);
for j = 1:length(VkClassRows)
    errArray(j,:,:) = errAllNormsMean{j,1}';
end

%%
% titlesNorm={'Norm 1 (and Norm Infinity)','Norm 2','Frobenius Norm'};
load('MatFiles\errKtest_K_2_2_40.mat');
ylabelNorm={'$e1$','$e_2$','$e_f$'};

for normI = 1:3
    figure;
    plot(VkClassRows,errArray(:,1,normI))
    hold on
    plot(VkClassRows,ones(length(VkClassRows),1)*errArray(1,2,normI))
    plot(VkClassRows,ones(length(VkClassRows),1)*errArray(1,3,normI))
    h=legend('$\widetilde{U}_{K}$ (Clustering)','$\widetilde{U}_{K}$ (Known Clustering)','${U}_{K}$ (PCA)')
    set(h,'Interpreter','Latex','Fontsize',12)
%     title(titlesNorm{normI});
    xlabel('K [Clusters]','Fontsize',15)
    ylabel(ylabelNorm{normI},'Fontsize',15)
    xlim([12 40])
end

