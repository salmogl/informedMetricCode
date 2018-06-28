clear all;
close all;
clc; 

dataFolderPathLoad  = '';
dataFolderPathSave  = '';
genesData           = xlsread(strcat(dataFolderPathLoad,'DCLungStudy_dChip-Processed_microarray_data'),5,'B2:DB22285');
genesData(:,35)     = []; % remove null
genesData(7206,:)   = []; % remove null

[~,subjectsStr]     = xlsread(strcat(dataFolderPathLoad,'DCLungStudy_dChip-Processed_microarray_data'),5,'B1:DB1');
subjectsStr(35)     = []; % remove null
[num,str]           = xlsread(strcat(dataFolderPathLoad,'DCLungStudy_Clinical_Covariates_with_Hgrade'),1);
str([37 104 105],:) = []; % remove null
num([36 103 104],:) = []; % remove null
month               = num(1:104,18);
vitalStat           = str(2:105,14);
strStageN           = str(2:105,21);
strStageT           = str(2:105,22);
subjectsStrProp     = str(2:105,5);

[M,N] = size(genesData);
% Sort subjects according to their DC_STUDY_ID - Genes
subjects        = zeros(N,1);
for ii = 1:N
    temp            = strsplit(subjectsStr{ii}(1:end-1),'133A_');
    temp            = strsplit(temp{2},'L');
    subjects(ii)    = str2double(temp{1});
end

[~,sortInd]             = sort(subjects);
dataStruct.genesData    = genesData(:,sortInd);

% Sort subjects according to their DC_STUDY_ID - Properties
subjects = zeros(N,1);
for ii = 1:N
    temp            = strsplit(subjectsStrProp{ii}(1:end-1),'133A_');
    temp            = strsplit(temp{2},'L');
    subjects(ii)    = str2double(temp{1});
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
    temp        = strsplit(strStageN{ii},'N');
    stageN      = str2double(temp{2}(1));
    temp        = strsplit(strStageT{ii},'T');
    stageT      = str2double(temp{2}(1));
    stage(ii)   = stageN*4 + stageT;
end

% VITAL_STATUS
dataStruct.vitalStat    = vitalStat01(sortIndProp);
% MONTHS_TO_LAST_CONTACT_OR_DEATH
dataStruct.month        = month(sortIndProp);
% PATHOLOGIC_STAGE
dataStruct.stage        = stage(sortIndProp);


save(strcat(dataFolderPathSave,'dataStruct_MSK'),'dataStruct')