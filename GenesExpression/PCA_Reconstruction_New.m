function [mDist,sumD,sd_k] = PCA_Reconstruction_New(mY,args)

rank               = args.rank;
nClust             = args.nClust;
knn                = args.knn;
replicates         = args.replicates;
informedFlag       = args.informedFlag;

[M,N]              = size(mY);
mDist              = zeros(N, N);
sumD               = 0;
sd_k               = 0;
if (knn < N)
   
    invCov = zeros(size(mY, 1), size(mY, 1), N);

    for ii = 1 : N
        ii
        mBurst_idx = knnsearch(mY', mY(:, ii)', 'K', knn);
    %     mBurst_idx_Clust = knnsearch(mY', mY(:, ii)', 'K', P_Clust);

    %     if ii < knn / 2
    %         mBurst_idx = 1 : knn;
    %     elseif ii >= size(mY, 2) - knn / 2 
    %         mBurst_idx = N - knn + 1: N;
    %     else
    %         mBurst_idx = ii - knn / 2 + 1 : ii + knn / 2;
    %     end
        mBurst_Y   = mY(:, mBurst_idx);
    %     mBurst_Y_Clust   = mY(:, mBurst_idx_Clust);

        %-- Center Data:
    %     mBurst_Y = bsxfun(@minus, mBurst_Y, mean(mBurst_Y, 2));

    %     nIter = 30;
    %     clust = zeros(M,nIter);
    %     sumd  = zeros(1,nIter);
    %     
    %     for jj = 1 : nIter
    %         [clust(:,jj),~,temp] = kmeans(mBurst_Y,nClust,'MaxIter',1e6,'OnlinePhase','on');
    % %           [clust(:,jj),~,temp]  = kmeans(mBurst_Y,rank+1,'Start',centroids(:,mBurst_idx),'MaxIter',1e6);
    %         h = histcounts(clust(:,jj),nClust);
    %         if(sum(find(h<2)) > 0)
    %             sumd(jj)             = inf;
    %         else
    %             sumd(jj)             = sum(temp);
    %         end
    %     end
    %     [~,ind] = min(sumd);
    %     clust   = clust(:,ind);

        if (informedFlag)
            [clust,~,D]   = kmeans(mBurst_Y,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',replicates);
        end
        %%%%%%%%%  Calculate Gap Statistic  %%%%%%%%%%%
%         if (informedFlag)
%             X             = mBurst_Y';
%             X             = X - repmat(mean(X,2),1,M);
%             nr            = histcounts(clust);
%             Wk            = 0.5*sum(D'./nr);
%             [~,~,Vb]      = svd(mBurst_Y');
%             X_            = X*Vb; 
%             maxLim        = repmat(max(X_,[],2),1,M);
%             minLim        = repmat(min(X_,[],2),1,M);
%             B             = 100;
%             Wkb           = zeros(1,B);
% 
%             for bb = 1:B
% 
%                 Z_            = minLim + (maxLim - minLim).*rand(knn,M);
%                 Z             = Z_*Vb';
%                 [~,~,Db]      = kmeans(Z,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',10);
%                 nrb           = histcounts(clust);
%                 Wkb(bb)       = 0.5*sum(Db'./nrb);
% 
%             end
% 
%             sumD = sumD + (mean(log(Wkb))-log(Wk))/N;
%             sd_k = std(log(Wkb));
%         end
        
        %sumD        = sumD + (log(knn*M/12) + 2*log(nClust)/knn- sum(D'./nr))/N;
        %sumD        = sum(D'./nr)/N + sumD;
        
    %     clust   = KmeansSvd(mBurst_Y',nClust,false,false,'PCA');

        if(knn > size(mBurst_Y,1))
            [~,S,V] = svd(cov(mBurst_Y'));
        else
            mBurst_Y = mBurst_Y - repmat(mean(mBurst_Y,2),1,knn);
            [~,S,V0] = svd(mBurst_Y'*mBurst_Y/(knn-1));
            V = mBurst_Y*V0*(S^-0.5)/sqrt(knn-1);
        end    
        
        if (informedFlag)
            invCov(:,:,ii) = calcInvCovInf(clust,S(1:rank,1:rank),rank,mBurst_Y); 
        else
            invCov(:,:,ii) = calcInvCovPca(V(:,1:rank),S(1:rank,1:rank),rank);
        end
%         nr          = histcounts(clust);
%         Wk          = 0.5*sum(D'./nr);
%         gap         = log(knn*M/12)-(2/knn)*log(nClust) - log(Wk)
    %     Ryy          = cov(mBurst_Y');
    %     tPCA(:,:,ii) = pinv(Ryy);
    end
    
else
    if (informedFlag)
        [clust,~,D] = kmeans(mY,nClust,'MaxIter',1e6,'OnlinePhase','on','Replicates',replicates);
    end
%     nr          = histcounts(clust);
%     sumD        = log(knn*M/12) + 2*log(nClust)/knn- sum(D'./nr);
%     sumD        = sum(D'./nr);
%         clust   = KmeansSvd(mY',nClust,false,false,'PCA');

    if(knn > M)
        [~,S,V] = svd(cov(mY'));
    else
        mY = mY - repmat(mean(mY,2),1,knn);
        [~,S,V0] = svd(mY'*mY/(knn-1));
        V = mY*V0*(S^-0.5)/sqrt(knn-1);
    end    
    
    if (informedFlag)
        tempInvCov = calcInvCovInf(clust,V(:,1:rank),rank,mY); 
    else
        tempInvCov = calcInvCovPca(V(:,1:rank),S(1:rank,1:rank),rank);
    end
            
    invCov = repmat(tempInvCov,[1 1 N]);
    
end

for ii = 1 : N
%     ii
    vY1 = mY(:,ii);
    for jj = ii + 1 : N
        vY2   = mY(:,jj);
        vDiff = vY1 - vY2;
        
            mDist(ii,jj) =  0.5 * vDiff' * ( invCov(:,:,ii) + invCov(:,:,jj) ) * vDiff;

    end
end

mDist = mDist + mDist';
mDist = sqrt(mDist); 
  
end
    
