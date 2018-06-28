clc;
clear all;
close all hidden;
close all;

set(0,'defaulttextinterpreter','latex')

nTrial          = 1;
nColsTest       = 100;
VkClassRows     = 18;
plotFlag        = true;
iterationType   = 'nSamples'; %'kClusters';'nSamples'

if (plotFlag)
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

switch iterationType
    
    case 'nSamples'
        nIters = length(nColsTest);
    case 'kClusters'
        nIters = length(VkClassRows);
        
end

% Data Generation
args.nRows              = 900;
args.kClassRows         = 18;
args.randFlag           = false;

% Affifnity parameters
paramsAff.eps           = 1;
paramsAff.metric        = 'euclidean';
paramsAff.thresh        = 0.05;

% Projected gradient parameters
paramsGrad.rank         = args.kClassRows;
paramsGrad.mu           = 3e-4;
paramsGrad.niter        = 200;
paramsGrad.stepOpt      = 'diminish';
paramsGrad.b            = 1;
paramsGrad.debugMode    = false;
paramsGrad.mu           = 1e-1;
paramsGrad.stepOpt      = 'constStepLength';
paramsGrad.threshDiff   = 0;

errAllNormsMean         = cell(nIters,2);
kClassRows              = args.kClassRows;
dataTypa                = 'noPerm';

for iter = 1:nIters;

    errAllNorms             = zeros(4,3);
    
    if(strcmp(iterationType,'nSamples'))
        args.nColsTest          = nColsTest(iter);
        paramsAff.knn           = nColsTest(iter);
    else
        args.nColsTest          = nColsTest(1);
        paramsAff.knn           = nColsTest(1);
    end
    
    for trial = 1:nTrial

        [data,bordersRows,covariance,columnsClass] = getData(args);

        data                = data - repmat(mean(data,2),1,size(data,2));
        [realPd,sReal,v]    = svd(covariance);
        [sizeRows,sizeCols] = size(data);

        switch(dataTypa)
            case 'noPerm'
                shuffleRows = 1:sizeRows;
            case 'randPerm'
                shuffleRows = randperm(sizeRows);
        end

        shuffleCols         = columnsClass;

        % Embedding computation using euclidean/cosine distance 
        aff = CalcInitAff(data,paramsAff);

        aff = stochastic(aff);
        eigsnum = 4;
        [eigvecs, eigvals] = eigs(aff, eigsnum);
        embedding = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

        % Embedding computation using Mahalanobis distance 
        [V,S,U]                      = svd(cov(data'));
        dataTestProjectedMahalanobis = (V(:,1:kClassRows-1)/sqrt(S(1:kClassRows-1,1:kClassRows-1)))'*data;
        aff                          = CalcInitAff(dataTestProjectedMahalanobis,paramsAff);
        aff                          = stochastic(aff);
        eigsnum                      = 4;
        [eigvecs, eigvals]           = eigs(aff, eigsnum);
        embeddingMahalanobis         = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

        % Embedding computation using Mahalanobis distance informed by clustering
        if(strcmp(iterationType,'kClusters'))
            curK        = VkClassRows(iter);
        else
            curK        = VkClassRows(1);
        end   
        [class,~,D]     = kmeans(data,curK,'OnlinePhase','on','Replicates',500);
        U_inf           = getPdSubProj(data,class,paramsGrad);
        S_hInv          =  diag(sum(S(1:kClassRows-1,1:kClassRows-1)).^-0.5);

        % project data on the informed principal directions
        dataTestProjectedClass  = (U_inf(:,1:kClassRows-1)*S_hInv)'*data;   
        
        aff                 = CalcInitAff(dataTestProjectedClass,paramsAff);
        aff                 = stochastic(aff);
        eigsnum             = 4;
        [eigvecs, eigvals]  = eigs(aff, eigsnum);
        embeddingKmeans     = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);

        [~,index]       = sort(shuffleRows);
        classRealRows   = zeros(sizeRows,1);
        
        for i=1:kClassRows
                classRealRows(index(bordersRows(i)+1:bordersRows(i+1))) = i;
        end

        U_infR          = getPdSubProj(data,classRealRows,paramsGrad);
        U_real          = realPd(:,1:kClassRows);
        U_pca           = V(:,1:kClassRows);
        S_Pca           = S(1:kClassRows,1:kClassRows);

        if (plotFlag)

            figure;
            imagesc(data), colormap jet, axis on

            covPCA   = U_pca*S_Pca*U_pca';
            covClass = U_inf*S_Pca*U_inf';

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
            h=legend('$\widetilde{\Sigma}_{i900}$','$\Sigma_{i900}$','$\widehat{\Sigma}_{i900}$');
            set(h,'Interpreter','latex')
            set(h,'Position',[0.2 0.7 0.2 0.2])
            set(h,'fontsize',15)
            ylim([-1 29])

            pointSize = 15;
            figure;
            scatter3(embeddingMahalanobis(:,1),embeddingMahalanobis(:,2),embeddingMahalanobis(:,3),pointSize,shuffleCols,'filed')
            xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
            ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
            zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
            view([-50 12])
            colormap winter

            figure;
            scatter3(embeddingKmeans(:,1),embeddingKmeans(:,2),embeddingKmeans(:,3),pointSize,shuffleCols,'filed')
            xlabel('$\phi_1$','Interpreter','latex','FontSize',15)
            ylabel('$\phi_2$','Interpreter','latex','FontSize',15)
            zlabel('$\phi_3$','Interpreter','latex','FontSize',15)
            view([-50 12])
            colormap winter


        end

        errNorm2Class           = norm(U_inf*U_inf'-U_real*U_real');
        errNorm2RealClass       = norm(U_infR*U_infR'-U_real*U_real');
        errNorm2Mahalanobis     = norm(U_pca*U_pca'-U_real*U_real');

        errNorm1Class           = norm(U_inf*U_inf'-U_real*U_real',1);
        errNorm1RealClass       = norm(U_infR*U_infR'-U_real*U_real',1);
        errNorm1Mahalanobis     = norm(U_pca*U_pca'-U_real*U_real',1);

        errNormInfClass         = norm(U_inf*U_inf'-U_real*U_real',inf);
        errNormInfRealClass     = norm(U_infR*U_infR'-U_real*U_real',inf);
        errNormInfMahalanobis   = norm(U_pca*U_pca'-U_real*U_real',inf);

        errNormFroClass         = norm(U_inf*U_inf'-U_real*U_real','fro');
        errNormFroRealClass     = norm(U_infR*U_infR'-U_real*U_real','fro');
        errNormFroMahalanobis   = norm(U_pca*U_pca'-U_real*U_real','fro');


        errAllNorms = errAllNorms +[errNorm1Class, errNorm1RealClass,errNorm1Mahalanobis;...
            errNorm2Class, errNorm2RealClass, errNorm2Mahalanobis;...
            errNormFroClass, errNormFroRealClass, errNormFroMahalanobis;...
            errNormInfClass, errNormInfRealClass, errNormInfMahalanobis;
            ];


    end

    errAllNormsMean{iter,1} = errAllNorms/nTrial;

    switch iterationType
    
    case 'nSamples'
        errAllNormsMean{iter,2} = strcat('Number of samples: ',num2str(nColsTest(iter)));
    case 'kClusters'
        errAllNormsMean{iter,2} = strcat('Number of clusters: ',num2str(VkClassRows(iter)));       
    end
    

end
set(0,'DefaultFigureVisible','on')

errArray = zeros(nIters,3,4);
for j = 1:nIters
    errArray(j,:,:) = errAllNormsMean{j,1}';
end

ylabelNorm={'$e1$','$e_2$','$e_f$'};

switch iterationType
    case 'nSamples'
        errAxis = nColsTest;
        labelAx = 'n [Samples]';
    case 'kClusters'
        errAxis = VkClassRows;
        labelAx = 'K [Clusters]';
end

for normI = 1:3
    figure;
    plot(errAxis,errArray(:,:,normI))
    h=legend('$\widetilde{U}_{K}$ (Clustering)','$\widetilde{U}_{K}$ (Known Clustering)','${U}_{K}$ (PCA)');
    set(h,'Interpreter','Latex','Fontsize',12)
    xlabel(labelAx,'Fontsize',15)
    ylabel(ylabelNorm{normI},'Fontsize',15)
end

%%

embeddingKmeansLabeled      = [embeddingKmeans shuffleCols];
embeddingMahalanobisLabeled = [embeddingMahalanobis shuffleCols];
[~, accuracyKmeans]         = trainClassifierR3(embeddingKmeansLabeled);
[~, accuracyMahalanobis]    = trainClassifierR3(embeddingMahalanobisLabeled);


