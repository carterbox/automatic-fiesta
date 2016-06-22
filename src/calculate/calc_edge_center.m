function data = calc_edge_center(data)
% Calculates the fractional distance from the center to the edge of the
% data for each point.

x = data.getPredictor('x');
y = data.getPredictor('y');

x0 = min(x(:));
x1 = max(x(:));

y0 = min(y(:));
y1 = max(y(:));

xleft = x - x0;

xright = x1 - x;

yleft = y - y0;
yright = y1 - y;

% Return the shortest manhattan distance to an edge
d = min(min(xleft,xright),min(yleft,yright));

radius = max(d(:));
r = (radius - d)/radius;

data = data.addPredictor('radius',r);
end
