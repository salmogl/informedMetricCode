clear all;
close all;
clc; 

addpath('..\mahalanobisEquivalence\3DQuest\3DQuest\Questionnaire')
addpath('kmplot')
load('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_MSK.mat');

data = dataStruct.genesData;
month = dataStruct.month;
vitalStat = dataStruct.vitalStat;

[val,indMaxStd] = sort(std(data,[],2),'descend');
data = data(indMaxStd(2:1000),:);
[M,N] = size(data);

% D = diag(sqrt(sum(data.^2,2)));
% data = D^-1*data;
% data = data'+ abs(min(data(:)));
% data = data/max(data(:));
% data = data*2-1;

figure;
imagesc(data)
title('Data')
ylabel('Genes')
xlabel('Subjects')

class = kmeans(data,15);
[val,indSort] = sort(class);

figure;
subplot(121)
imagesc(data(indSort,:))
title('Data')
ylabel('Clustered Genes')
xlabel('Subjects')
subplot(122)
imagesc(class(indSort,:))
title('Clusters')

figure;
imagesc(cov(data(indSort,:)'))
title('Subjects Covariance')
% view([-54,-51])
data = data - repmat(mean(data,2),1,N);
[U,S,V] = svd(cov(data'));

classNum = length(unique(class));
Classifier = zeros(M, classNum);
for ii = 1:classNum
    Classifier(class == ii,ii) = 1;
end

% normalization according to group size
ClassifierNorm = Classifier./(repmat(sqrt(sum(Classifier)),[M 1]));
Q = [ones(M,1)/M U(:,1:classNum-1)];
T = ((ClassifierNorm'*ClassifierNorm)^-1)*ClassifierNorm'*Q;
T = T(:,2:end)/sqrt(S(1:classNum-1,1:classNum-1));
pcClass = ClassifierNorm*T;

pcClassNorm = pcClass./repmat(sqrt(sum(pcClass.^2,1)),[M 1]);

dataProjClust = (pcClassNorm/sqrt(S(1:classNum-1,1:classNum-1)))'*data;
dataProjPca = (V(:,1:classNum-1)/sqrt(S(1:classNum-1,1:classNum-1)))'*data;
% dataProjClust = pcClassNorm'*data;
% dataProjPca = V(:,1:classNum-1)'*data;

paramsAff.eps = 5;
paramsAff.knn = 20;
paramsAff.metric = 'euclidean'; % 'cosine_similarity'/'euclidean'
paramsAff.thresh = 0.05;
aff = CalcInitAff(dataProjPca,paramsAff);
aff = stochastic(aff);

eigsnum = 10;
[eigvecs, eigvals] = eigs(aff, eigsnum);
embeddingPca= eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

color = jet(max(month)+1);
colorsVital = color(month,:);
colorsVital(vitalStat == 1,:) = repmat(color(max(month)+1,:),length(find(vitalStat)), 1);
colorsStage = dataStruct.stage; %floor(dataStruct.stage/4);

figure;
subplot(221)
plot3dData(embeddingPca(:,1:3)',colorsStage,'Embedding Colored By Stage (PCA)')
subplot(222)
plot3dData(embeddingPca(:,1:3)',colorsVital,'Embedding Colored By Month (PCA)')
subplot(223)
plot3dData(dataProjPca(1:3,:),colorsStage,'First 3 PC Colored By Stage (PCA)')
subplot(224)
plot3dData(dataProjPca(1:3,:),colorsVital,'First 3 PC Colored By Month (PCA)')

aff = CalcInitAff(dataProjClust,paramsAff);
aff = stochastic(aff);
[eigvecs, eigvals] = eigs(aff, eigsnum);
embeddingClust= eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

figure;
subplot(221)
plot3dData(embeddingClust(:,1:3)',colorsStage,'Embedding Colored By Stage (Clustering)')
subplot(222)
plot3dData(embeddingClust(:,1:3)',colorsVital,'Embedding Colored By Month (Clustering)')
subplot(223)
plot3dData(dataProjClust(1:3,:),colorsStage,'First 3 PC Colored By Stage (Clustering)')
subplot(224)
plot3dData(dataProjClust(1:3,:),colorsVital,'First 3 PC Colored By Month (Clustering)')



% groupsPca = kmeans(embeddingPca,2);
% groupsClust = kmeans(embeddingClust,2);
groupsPca = KmeansSvd(embeddingPca',2,true,true,'PCA');
groupsClust = KmeansSvd(embeddingClust',2,true,true,'Clustering');



[vHistPca0, tAxisPca0] = getVitalHist(vitalStat(groupsPca == 1),month(groupsPca == 1));
[vHistPca1, tAxisPca1] = getVitalHist(vitalStat(groupsPca ~= 1),month(groupsPca ~= 1));
% [vHistPca2, tAxisPca2] = getVitalHist(vitalStat(groupsPca == 3),month(groupsPca == 3));

figure;
[~,~,a1,b1] = kmplot([ month(groupsPca == 1) , vitalStat(groupsPca == 1) ]);
[~,~,a2,b2] = kmplot([ month(groupsPca ~= 1) , vitalStat(groupsPca ~= 1) ]);

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

subplot(122)
stairs(a1,b1,'r')
hold on
stairs(a2,b2,'b')
xlabel('Month')
ylabel('Proportion Alive')
xlim([1 60])
ylim([0 1])
title('PCA (kmplot)')

[vHistClust0, tAxisClust0] = getVitalHist(vitalStat(groupsClust == 1),month(groupsClust == 1));
[vHistClust1, tAxisClust1] = getVitalHist(vitalStat(groupsClust ~= 1),month(groupsClust ~= 1));
% [vHistClust2, tAxisClust2] = getVitalHist(vitalStat(groupsClust == 3),month(groupsClust == 3));

figure;
[~,~,a1,b1] = kmplot([ month(groupsClust == 1) , vitalStat(groupsClust == 1) ]);
[~,~,a2,b2] = kmplot([ month(groupsClust ~= 1) , vitalStat(groupsClust ~= 1) ]);

figure;
subplot(121)
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

subplot(122)
stairs(a1,b1,'r')
hold on
stairs(a2,b2,'b')
xlabel('Month')
ylabel('Proportion Alive')
xlim([1 60])
ylim([0 1])

title('Clustering (kmplot)')





% figure;
% scatter3(VClass(:,1),VClass(:,2),VClass(:,3),pointSize,vitalStat+1,'filed')



