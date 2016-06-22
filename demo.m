addpath(genpath('./src'));
addpath('./models');
addpath(genpath('/home/chingd/Documents/MATLAB/wpart'));

%% A little bit about the files.
JSON_dir = '/home/chingd/Documents/MATLAB/strains_scripts/Cooktime/JSON';
MAT_dir = '/home/chingd/Documents/MATLAB/strains_scripts/Cooktime/MAT';

%% Collect all of the statistical data.

for file = dir(JSON_dir)'
    if ~file.isdir     
        constructStatsData([file.folder '/' file.name], MAT_dir);
    end
end

%% Collect all of the statistical data in to one object.

alldata = cell(5,1); i = 1;
for file = dir(MAT_dir)'
    if ~file.isdir
        afilename = [file.folder '/' file.name];
        load(afilename)
        
        disp(i);
        alldata{i} = data;
        i = i + 1;
    end
end

data = statsdata();
data = data.combine(alldata);

clear alldata

%% Fit a model to the data.

[response, rheader] = data.getResponse('Em'); % Choose response varible
[predictors, pheaders] = data.getPredictors(0); % Provide all predictors
headers = [rheaders pheaders]; % concatenate the data labels

betas = [0, 0]; % provide initial guesses for the parameters
functiontemplate = 'y ~ PretendModel(b1, b2, 42, x1, x2, x3)';

% fit the model
M = fitnlm(predictors, response, functiontemplate, betas,...
           'VarNames', headers, 'CoefficientNames', bheaders);