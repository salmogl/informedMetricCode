clc;
clear all;
close all;


dataSet  = 'CANDF';
loadPath = addPathLoad(dataSet);
load(loadPath);

data            = dataStruct.genesData;
[~,indMaxStd]   = sort(std(data,[],2),'descend');
data            = data(indMaxStd(2:201),:);
[M,N]           = size(data);
nClust          = 2:15;%12:2:24;
sumD            = zeros(1,length(nClust));
knn             = 34;
parfor ii = 1:length(nClust)
ii
for jj =1:N
jj
mBurst_idx      = knnsearch(data', data(:, jj)', 'K', knn);
mBurst          = data(:,mBurst_idx);
[clust,~,D]     = kmeans(mBurst,nClust(ii),'MaxIter',1e6,'OnlinePhase','on','Replicates',1e2);
X               = mBurst';
X               = X - repmat(mean(X,2),1,M);
nr              = histcounts(clust);
Wk              = 0.5*sum(D'./nr);
[~,~,Vb]        = svd(mBurst');
X_              = X*Vb; 
maxLim          = repmat(max(X_,[],2),1,M);
minLim          = repmat(min(X_,[],2),1,M);
B               = 1;
Wkb             = zeros(1,B);

    for bb = 1:B

        Z_            = minLim + (maxLim - minLim).*rand(knn,M);
        Z             = Z_*Vb';
        [~,~,Db]      = kmeans(Z,nClust(ii),'MaxIter',1e6,'OnlinePhase','on','Replicates',50);
        nrb           = histcounts(clust);
        Wkb(bb)       = 0.5*sum(Db'./nrb);

    end

sumD(ii) = sumD(ii) + (mean(log(Wkb))-log(Wk))/N;
% sd_k     = std(log(Wkb));
end
end

figure;
plot(nClust,sumD);

%%
[val,indSort] = sort(clust);
figure;
HMObject = HeatMap(data(indSort,:));