function [f] = errorbar_plot(data, predictor, responses, f0)
% Plots and line plot with error bars for each of the responses as a
% function of the predictor on one plot plot.
% ------------------------------------------------------------------------ 

[X, Xheaders] = data.getPredictor(predictor);
[Y, Yheaders] = data.getResponse(responses);

groups_of_x = unique(X);

means_of_y = zeros(numel(groups_of_x), size(Y,2));
error_of_y = means_of_y;

% calculate means and errorbars
for i = 1:numel(groups_of_x)
    means_of_y(i,:) = mean(Y(X == groups_of_x(i),:), 1);
    error_of_y(i,:) =  std(Y(X == groups_of_x(i),:), 1);
end

% plot the errorbars and means
if nargin > 3
    f = figure(f0); hold on;
else
    f = figure(); hold on;
end

for y = 1:size(Y,2)
    errorbar(groups_of_x, means_of_y(:,y), error_of_y(:,y));
end
xlabel(Xheaders)
legend(Yheaders);
hold off;
end