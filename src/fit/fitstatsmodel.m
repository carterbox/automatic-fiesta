function [model, report] = fitstatsmodel(data, model, betas, predictors, response, randomsample)
%FITSTATSMODEL fits the coefficients fot the given model, using predictors
% (strings) and response(string).
% betas: an Nx3 cell with headers, initalvalue, units

[pred,predh,predu] = data.getPredictor(predictors);
[resp,resph,respu] = data.getResponse(response);
beta = [betas{:,2}];
betah = {betas{:,1}};
betau = {betas{:,3}};

% take a random sample
if randomsample
    k = floor(randomsample*data.ocount);
    i = randperm(data.ocount,k);
    pred = pred(i,:);
    resp = resp(i,:);
end

model = fitnlm(pred, resp, model, beta,...
               'VarNames',[predh, resph],...
               'CoefficientNames', betah);
beta = model.Coefficients.Estimate;
report = 0;
end

