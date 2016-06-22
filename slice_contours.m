function f2 = slice_contours(data, step, predictors, response)
% creates a series of 3D scatter plots at each step. Where the the z axis
% of the scatter plot is the response and the x,y axes are the predictors.

[S, Sheader] = data.getPredictor(step);
[XY, XYheaders] = data.getPredictor(predictors);
[Z, Zheader] = data.getResponse(response);

for s in S
   
    
end


end

