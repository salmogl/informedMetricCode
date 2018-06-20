clc;
clear all;
close all;
set(0,'defaulttextinterpreter','latex');


% Generate data in the original space X (in R^2)
N         = 500; 
M         = 1200;       % #Fetures in high dimensional space
x1        = rand(N,1);
x2        = rand(N,1);
X1        = repmat(x1',[M 1]);
X2        = repmat(x2',[M 1]);
dataOrig  = [x1' ; x2'];

% parameters for the linear transformation
sigma     = 0.01*[1 0;0 1];
mu1       = [1/3 1];       
mu2       = [1/3 -1];         
mu3       = 3*[-1/3 -1];  
transform = [mvnrnd(mu1,sigma,M/3) ; mvnrnd(mu2,sigma,M/3) ; mvnrnd(mu3,sigma,M/3)];

% tranform to high dimensional space
data            = transform*dataOrig;
data            = data - repmat(mean(data,2),1,size(data,2));

% compute diffusion maps for the samples in X-space
args.eps        = 1;
distOrig        = pdist(dataOrig');
embeddingOrig   = calcEmbedding(distOrig,args);

%%
% Rows clustering
plotFlagKmeans = false;
useFit         = true;
title_         = ' Rows';
threshold      = 0.05;
k              = 3;
class          = kmeans(data,k);

[U,S,~]        = svd(cov(data'));
[~,indSort]    = sort(class);

% organize the cluster solution in the matrix H
classNum   = length(unique(class));
Classifier = zeros(size(class,1), classNum);
for i=1:classNum
    Classifier(class == i,i) = 1;
end
% normalization according to group size
H = Classifier./(repmat(sqrt(sum(Classifier)),[size(class,1) 1]));

% compute the informed principal directions 
paramsGrad.rank      = 2;
paramsGrad.mu        = 3e-4;
paramsGrad.niter     = 100;
paramsGrad.stepOpt   = 'diminish';
paramsGrad.b         = 1;
paramsGrad.mu        = 1e-1;
paramsGrad.stepOpt   = 'constStepLength';
paramsGrad.debugMode = false;
paramsGrad.threshDiff= 0;
U_c                  = getPdSubProj(data,class,paramsGrad);


S_hInv              =  diag(sum(S(1:2,1:2)).^-0.5);

% project data on standart pc and informed pc
dataProjectedClass  = (U_c*S_hInv)'*data;
dataProjected       = (U(:,1:2)*S_hInv)'*data;


%%
% compute distances and embeddings
mahabLikeDistClass  = pdist(dataProjectedClass');
embeddingClass      = calcEmbedding(mahabLikeDistClass,args);
embeddingClassDist  = pdist(embeddingClass);

%%
% correlation between the original dustances and the mahalanobis distances
distData            = pdist(data');
M                   = corrcoef([distOrig' mahabLikeDistClass']);
classCor            = M(2,1);
M                   = corrcoef([distOrig' distData']);
eucCor              = M(2,1);

mahabLikeDistClass  = mahabLikeDistClass/max(mahabLikeDistClass);
distOrig            = distOrig/max(distOrig);
distData            = distData/max(distData);
%% Figures

figure;
imagesc(data)
colormap(autumn(5))

pointSize           = 12;
figure;
scatter(embeddingClass(:,1),embeddingClass(:,2),pointSize,x1,'filed')
xlabel('$\phi_1$','FontSize',15)
ylabel('$\phi_2$','FontSize',15)

figure;
scatter(x1,x2,pointSize,cos(x1),'filed')
xlabel('$x_1$','Interpreter','latex','FontSize',15)
ylabel('$x_2$','Interpreter','latex','FontSize',15)

figure;
scatter(distOrig,mahabLikeDistClass,pointSize,'filled')
xlabel('$\|x_1-x_2\|_2$','FontSize',15)
ylabel('$\tilde{d}_{M}(c_1,c_2)$','Interpreter','latex','FontSize',15)
t=text(0.2,0.8,strcat('Correlation = ',sprintf('%.2f',classCor)));
t.FontSize = 15;

figure;
scatter(distOrig,distData,pointSize,'filled')
xlabel('$\|x_1-x_2\|_2$','FontSize',15)
ylabel('$\|c_1-c_2\|_2$','FontSize',15)
t=text(0.2,0.8,strcat('Correlation = ',sprintf('%.2f',eucCor)));
t.FontSize = 15;