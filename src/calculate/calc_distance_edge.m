function data = calc_distance_edge(data,x0,y0,width,height)
% Calculates the shortest manhattan distance from all the points to the edge.
% Edge as defined by the box whose min corner is x0 y0 and width and height.
% All the points in data should be insided the box, else there will be negative
% distances.

[x, ~, unit] = data.getPredictor('x');
y = data.getPredictor('y');

y1 = y0 + height;
x1 = x0 + width;

ym = y - mean([y0,y1]);
xm = x - mean([x0,x1]);

xleft = x - x0;
xright = x1 - x;

yleft = y - y0;
yright = y1 - y;

% Return the shortest manhattan distance to an edge
d = min(min(xleft,xright),min(yleft,yright));

data = data.addPredictor('edge',d,unit{1});
data = data.addPredictor('center',min(xm,ym),unit{1});
data = data.addPredictor('xcenter',xm,unit{1});
data = data.addPredictor('ycenter',ym,unit{1});
end
