function [model, report] = fitlmstatsmodel(data, type, predictors, response, randomsample)
%FITSTATSMODEL fits the coefficients fot the given model, using predictors
% (strings) and response(string).
% betas: an Nx3 cell with headers, initalvalue, units

[pred,predh,predu] = data.getPredictor(predictors);
[resp,resph,respu] = data.getResponse(response);

% take a random sample
if randomsample
    k = floor(randomsample*data.ocount);
    i = randperm(data.ocount,k);
    pred = pred(i,:);
    resp = resp(i,:);
end

model = stepwiselm(pred, resp, type,...
               'VarNames',[predh, resph],'PEnter', 1e-10, 'PRemove', 1e-2,...
               'CategoricalVars', [5,6]);
report = 0;
end

