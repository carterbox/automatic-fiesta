function data = calc_distance_center(data, center)
% Calculates distance from center predictor

x = data.getPredictor('x');
y = data.getPredictor('y');

r = sqrt((x-center(1)).^2 + (y-center(2)).^2);

data = data.addPredictor('radius',r);

end
