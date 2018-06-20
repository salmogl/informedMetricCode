function loadPath = addPathLoad(dataSet)

dir = strsplit(pwd,{'/','\'});

if (strcmp(dir{2},'home'))

	addpath('../../localMahalanobis');
	addpath('../../3DQuest/3DQuest/Questionnaire');
	addpath('kmplot');
	addpath('logrank');
	addpath('../');
    addpath('C:\Users\salmogl\Google Drive\Master\MATLAB\Matlab_Utils\tSNE\')

    
    loadPath = strcat('../../geneExpression/data_',dataSet,'.mat');

else

	addpath('..\..\..\localMahalanobis');
	addpath('..\..\..\mahalanobisEquivalence\3DQuest\3DQuest\Questionnaire');
	addpath('kmplot');
	addpath('logrank');
	addpath('..\');
    addpath('C:\Users\salmogl\Google Drive\Master\MATLAB\Matlab_Utils\tSNE\')

    loadPath = strcat('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_',dataSet,'.mat');

end

	