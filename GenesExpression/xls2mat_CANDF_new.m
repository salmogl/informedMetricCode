clear all;
close all;
clc; 


genesData = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data',4,'B2:CE22284');
% genesData =
% xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\R_new\my_combat_edata',1,'B2:CE501');
% COMBAT data


[num,subjectsStr] = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_dChip-Processed_microarray_data',4,'B1:CE1');
[num,str] = xlsread('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\DCLungStudy_Clinical_Covariates_with_Hgrade',1);
month = num(105:186,18);
vitalStat = str(106:187,14);
strStageN = str(106:187,21);
strStageT = str(106:187,22);
subjectsStrProp = str(106:187,5);


[M,N] = size(genesData);

% Sort subjects according to their DC_STUDY_ID - Genes
subjects = zeros(N,1);
for ii = 1:N
    temp = strsplit(subjectsStr{ii},'CL');
    temp = strsplit(temp{2},'A');
    subjects(ii) = str2double(temp{1});
end
[val,sortInd] = sort(subjects);
dataStruct.genesData = genesData(:,sortInd);

% Sort subjects according to their DC_STUDY_ID - Properties
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

% Get PATHOLOGIC_N_STAGE and PATHOLOGIC_T_STAGE
stage = zeros(N,1);
for ii = 1:N
    temp   = strsplit(strStageN{ii},'N');
    stageN = str2double(temp{2}(1));
    temp   = strsplit(strStageT{ii},'T');
    stageT = str2double(temp{2}(1));
    stage(ii) = stageN*4 + stageT;
end

% VITAL_STATUS
dataStruct.vitalStat = vitalStat01(sortIndProp);
% MONTHS_TO_LAST_CONTACT_OR_DEATH
dataStruct.month = month(sortIndProp);
% PATHOLOGIC_STAGE
dataStruct.stage = stage(sortIndProp);


save('C:\Users\salmogl\Documents\Data\SurvivalPredictGeo\Matfiles\data_CANDF.mat','dataStruct')