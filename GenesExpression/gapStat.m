function [sumD,sd_k] = gapStat(mY,args)

nClust             = args.nClust;
knn                = args.knn;
replicates         = args.replicates;

[M,N]              = size(mY);
sumD               = 0;
sd_k               = 0;

  
for ii = 1 : N

    mBurst_idx      = knnsearch(mY', mY(:, ii)', 'K', knn);
    mBurst_Y        = mY(:, mBurst_idx);
    [clust,~,D]     = kmeans(mBurst_Y,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',replicates);

    %%%%%%%%%  Calculate Gap Statistic  %%%%%%%%%%%
    X               = mBurst_Y';
    X               = X - repmat(mean(X,2),1,M);
    nr              = histcounts(clust);
    Wk              = 0.5*sum(D'./nr);
    [~,~,Vb]        = svd(mBurst_Y');
    X_              = X*Vb; 
    maxLim          = repmat(max(X_,[],2),1,M);
    minLim          = repmat(min(X_,[],2),1,M);
    B               = 100;
    Wkb             = zeros(1,B);

    for bb = 1:B

        Z_            = minLim + (maxLim - minLim).*rand(knn,M);
        Z             = Z_*Vb';
        [~,~,Db]      = kmeans(Z,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',10);
        nrb           = histcounts(clust);
        Wkb(bb)       = 0.5*sum(Db'./nrb);

    end

    sumD = sumD + (mean(log(Wkb))-log(Wk))/N;
    sd_k = std(log(Wkb));

end
     
    
