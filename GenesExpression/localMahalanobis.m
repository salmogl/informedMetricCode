function [mDistPca,mDistClust] = localMahalanobis(mY,args)

rank               = args.rank;
nClust             = args.nClust;
knn                = args.knn;
replicates         = args.replicates;

[M,N]              = size(mY);
mDistPca           = zeros(N, N);
mDistClust         = zeros(N, N);


if (knn < N)

    invCovPca   = zeros(size(mY, 1), size(mY, 1), N);
    invCovClust = zeros(size(mY, 1), size(mY, 1), N);

    for ii = 1 : N
        
        mBurst_idx = knnsearch(mY', mY(:, ii)', 'K', knn);
        mBurst_Y   = mY(:, mBurst_idx);
        clust      = kmeans(mBurst_Y,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',replicates);

        if(knn > size(mBurst_Y,1))
            [~,S,V] = svd(cov(mBurst_Y'));
        else
            mBurst_Y = mBurst_Y - repmat(mean(mBurst_Y,2),1,knn);
            [~,S,V0] = svd(mBurst_Y'*mBurst_Y/(knn-1));
            V = mBurst_Y*V0*(S^-0.5)/sqrt(knn-1);
        end
        
        [invCovClust(:,:,ii),invCovPca(:,:,ii)] = calcInvCov(clust,V(:,1:rank),S(1:rank,1:rank),rank,mBurst_Y);

    end
    
else
    
    clust   = kmeans(mY,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',replicates);

    if(knn > M)
        [~,S,V] = svd(cov(mY'));
    else
        mY = mY - repmat(mean(mY,2),1,knn);
        [~,S,V0] = svd(mY'*mY/(knn-1));
        V = mY*V0*(S^-0.5)/sqrt(knn-1);
    end
    
    [tempClust,tempPca] = calcInvCov(clust,V(:,1:rank),S(1:rank,1:rank),rank,mY);
    invCovClust         = repmat(tempClust,[1 1 N]);
    invCovPca           = repmat(tempPca,[1 1 N]);
    
end

for ii = 1 : N

    vY1 = mY(:,ii);
    
    for jj = ii + 1 : N
        vY2   = mY(:,jj);
        vDiff = vY1 - vY2;

        mDistPca(ii,jj)   = ...
                           0.5 * vDiff' * ( invCovPca(:,:,ii) + invCovPca(:,:,jj) ) * vDiff;
        mDistClust(ii,jj) = ...
                           0.5 * vDiff' * ( invCovClust(:,:,ii) + invCovClust(:,:,jj) ) * vDiff;
    end
end

mDistPca = mDistPca + mDistPca';
mDistPca = sqrt(mDistPca);

mDistClust = mDistClust + mDistClust';
mDistClust = sqrt(mDistClust); 
  
end
    
