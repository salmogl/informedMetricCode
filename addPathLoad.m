function loadPath = addPathLoad()

dir = strsplit(pwd,{'/','\'});

if (strcmp(dir{2},'home'))

	addpath('../../localMahalanobis')
	addpath('../../3DQuest/3DQuest/Questionnaire')
	addpath('kmplot')
	addpath('logrank')
	addpath('../')

	loadPath = '../../geneExpression/data_CANDF.mat';

else

	addpath('..\..\localMahalanobis')
	addpath('..\..\mahalanobisEquivalence\3DQuest\3DQuest\Questionnaire')
	addpath('kmplot')
	addpath('logrank')
	addpath('..\')

	loadPath = 'C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_CANDF.mat';
	
end

	