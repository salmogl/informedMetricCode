function U1 = getPdSubProj(X,clust,params)

rank    = params.rank;
mu      = params.mu;
niter   = params.niter;
b       = params.b;

sigma   = cov(X');
[U,~,~] = svd(sigma);
U       = U(:,1:rank);
m       = size(U,1);

% organize the cluster solution in the matrix H
classNum   = length(unique(clust));
Classifier = zeros(size(clust,1), classNum);
for i=1:classNum
    Classifier(clust == i,i) = 1;
end
% normalization according to group size
H = Classifier./(repmat(sqrt(sum(Classifier)),[size(clust,1) 1]));

U1       = U;

% U0       = rand(m,rank);
diff     = zeros(niter,1);


for ii = 1:niter

    U0       = U1;
    Grad     = (eye(m)-U1*U1')*sigma*U1 + sigma*(eye(m)-U1*U1')*U1;
%     Grad     = (eye(m)-U0*U0')*cov(X')*U0;
%     U0       = U0 + mu*Grad;%/norm(Grad,2);
    switch params.stepOpt
        
        case 'diminish'
            
            U1       = U1 + mu*Grad/sqrt(ii);
            
        case 'constStepLength'
            
            U1       = U1 + mu*Grad/norm(Grad,2);
            
        case 'constStepSize'
            
            U1       = U1 + mu*Grad;
            
        case 'sqareSummable'
            
            U1       = U1 + mu*Grad/(b+ii);
            
    end


    U1       = (H*H')*U1;
    
    diff(ii)  = norm(U0-U1,'fro');
    
    if (diff(ii) < params.threshDiff)
        break
    end
    
%     % debug
%     if (params.debugMode)
%         diff(ii)  = norm(U0-U1,'fro');%norm(U*U' - U0*U0'); 
%     end
    
    
end

% norm2(X-(U+a*(X-U*U'*X)*X'*U)*(U+a*(X-U*U'*X)*X'*U)'*X)^2
U1  = U1./repmat(sqrt(sum(U1.^2,1)),m,1);

% debug
if (params.debugMode)
    figure;
    plot(diff);
    xlabel('Iterations')
    ylabel('Error')
    ylim([0 1e-2]);
end