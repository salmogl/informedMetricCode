function [dataTest,bordersRows,covariance,columnsClassTest] = getData(args)


nRows       = args.nRows;
nColsTest   = args.nColsTest;
kClassRows  = args.kClassRows;
randFlag    = args.randFlag;

if (randFlag)
    bordersRows = 1 + randperm(nRows-2);
    bordersRows = bordersRows(1:kClassRows-1);
    bordersRows = [ 1 sort(bordersRows) nRows];
else
    bordersRows = [0 50 120 190 220 285 320 355 400 470 500 545 603 642 700 735 800 850 nRows];
end

dataTempTest    = zeros(nRows,nColsTest);
kClassRows      = length(bordersRows)-1;

x               = linspace(0.2,1,kClassRows);
muF             = 5 + 5*x.^4;
muM             = 5 - 5*x.^4;
sigma           = 1;
ro              = 0.6*sigma;
muFfull         = zeros(nRows,1);
muMfull         = zeros(nRows,1);
realCovariance  = zeros(nRows);

mOrFtest        = randi(2,[nColsTest 1]);
nMaleTest       = length(find(mOrFtest == 1));
nFemaleTest     = nColsTest - nMaleTest;

for classRowsInd = 1:kClassRows

rowsClassIndex  = bordersRows(classRowsInd)+1:bordersRows(classRowsInd+1);
rowsClassSize   = length(rowsClassIndex);
cov0            = sigma*eye(rowsClassSize);
cov0(cov0 == 0) = ro;

dataTempTest(rowsClassIndex,mOrFtest == 1) = ...
   mvnrnd(muF(classRowsInd)*ones(rowsClassSize,1),cov0,nMaleTest)'; 
dataTempTest(rowsClassIndex,mOrFtest == 2) = ...
   mvnrnd(muM(classRowsInd)*ones(rowsClassSize,1),cov0,nFemaleTest)'; 

muFfull(rowsClassIndex) = muF(classRowsInd);
muMfull(rowsClassIndex) = muM(classRowsInd);
realCovariance(rowsClassIndex,rowsClassIndex) = cov0;

end

covariance                  = realCovariance +(muFfull-muMfull)*(muFfull-muMfull)'/4;
dataTest                    = zeros(size(dataTempTest));
dataTest(:,1:nMaleTest)     = dataTempTest(:,mOrFtest == 1);
dataTest(:,nMaleTest+1:end) = dataTempTest(:,mOrFtest == 2);
columnsClassTest            = sort(mOrFtest);
