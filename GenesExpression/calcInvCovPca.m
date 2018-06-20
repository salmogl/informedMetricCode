function invCov = calcInvCovPca(pc,eignVals,rank)
   

sInv         = diag(1./sum(eignVals(1:rank,1:rank)));

invCov       = pc(:,1:rank)*sInv*pc(:,1:rank)';
