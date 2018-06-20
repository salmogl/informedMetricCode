clear all;
close all;
clc; 

% data=xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data','B2:CF22285');
% data=xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data',2,'B3:GL22285');

genesData = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data',4,'B2:CE22284');
[num,subjectsStrData] = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data',4,'B1:CE1');
[num,vitalStat] = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_Clinical_Covariates_with_Hgrade',1,'N106:N187');
[num,subjectsStrProp] = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_Clinical_Covariates_with_Hgrade',1,'E106:E187');
month = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_Clinical_Covariates_with_Hgrade',1,'R106:R187');


[M,N] = size(genesData);

% Sort subjects according to their DC_STUDY_ID
subjects = zeros(N,1);
for ii = 1:N
    temp = strsplit(subjectsStrData{ii},'CL');
    temp = strsplit(temp{2},'A');
    subjects(ii) = str2double(temp{1});
end
[val,sortInd] = sort(subjects);
dataStruct.genesData = genesData(:,sortInd);

% Sort subjects according to their DC_STUDY_ID
subjects = zeros(N,1);
for ii = 1:N
    temp = strsplit(subjectsStrProp{ii},'CL');
    temp = strsplit(temp{2},'A');
    subjects(ii) = str2double(temp{1});
end
[val,sortIndProp] = sort(subjects);

% Get vital status
vitalStat01 = zeros(N,1);
for ii = 1:N 
    if(strcmp('Alive',vitalStat{ii}))
        vitalStat01(ii) = 1;
    end
end

% MONTHS_TO_LAST_CONTACT_OR_DEATH
dataStruct.month = month(sortIndProp);
% VITAL_STATUS
dataStruct.vitalStat = vitalStat01(sortIndProp);

save('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_CANDF.mat','dataStruct')