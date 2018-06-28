function [invCovClass,invCovPseudo] = calcInvCov(class,pc,eignVals,rank,data)

sizeRows = size(class,1);
params.rank         = rank;
params.mu           = 1e-1;
params.b            = 1;
params.niter        = 500;
params.stepOpt      = 'constStepLength';
params.debugMode    = false;
params.threshDiff   = 5e-8; %5e-5

pcClass             = getPdSubProj(data,class,params);
pcClass             = pcClass./repmat(sqrt(sum(pcClass.^2,1)),sizeRows,1);
sInv                = diag(1./sum(eignVals(1:rank,1:rank)));
invCovClass         = pcClass(:,1:rank)*sInv*pcClass(:,1:rank)';
invCovPseudo        = pc(:,1:rank)*sInv*pc(:,1:rank)';
