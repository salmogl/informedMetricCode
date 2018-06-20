function [invCovClass,invCovPseudo] = calcInvCovAriel(class,pc,eignVals,rank,data)

% sizeRows = size(class,1);
    
classNum = length(unique(class));
Classifier = zeros(size(class,1), classNum);
for i=1:classNum
    Classifier(class == i,i) = 1;
end

% normalization according to group size

ClassifierNorm = Classifier./(repmat(sqrt(sum(Classifier)),[size(class,1) 1]));
% % Q = [ones(sizeRows,1)/sqrt(sizeRows) pc(:,1:classNum-1)];
% Q = [ones(sizeRows,1)/sqrt(sizeRows) pc(:,1:rank)];
% T = ((ClassifierNorm'*ClassifierNorm)^-1)*ClassifierNorm'*Q;
% % T = T(:,2:end)/sqrt(eignVals(1:classNum-1,1:classNum-1));
% % T = T(:,2:end);
% T   = T(:,2:end);
% pcClass = ClassifierNorm*T;
% normalize pc
% pcClass = pcClass./repmat(sqrt(sum(pcClass.^2,1)),sizeRows,1);
params.rank      = rank;
% params.mu        = 5e-9;
params.mu        = 1e-1;
params.b         = 1;
params.niter     = 500;
params.stepOpt   = 'constStepLength';
% pcClass     = getPdSubProj(data,class,params);
[~,S_infR,V_infR] = svd(cov((ClassifierNorm*ClassifierNorm'*data)'));
pcClass = V_infR(:,1:rank);
% S_infR  = S_infR(1:rank,1:rank);
% classNum = 5;

% sInv = diag(1./sum(eignVals(1:classNum-1,1:classNum-1)));
sInv  = diag(1./sum(eignVals(1:rank,1:rank)));
% sInvR = diag(1./sum(S_infR(1:rank,1:rank)));

invCovClass = ...
      pcClass(:,1:rank)*sInv*pcClass(:,1:rank)';
%     pcClass(:,1:classNum-1)*sInv*pcClass(:,1:classNum-1)';

invCovPseudo = ...
      pc(:,1:rank)*sInv*pc(:,1:rank)';
%     pc(:,1:classNum-1)*sInv*pc(:,1:classNum-1)';