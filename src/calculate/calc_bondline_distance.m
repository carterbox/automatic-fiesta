function data = calc_bondline_distance(data, bondline_file)
% Calculate euclidian distance from each point to the bondline.
%
% Fits a bond plane from points provided in a file. Then from this best
% fit plane, bondline distances are calculated for each of the points in
% the statsdata. Negative distances are on the min side of the bondline
% and positive distances are on the max side of the bondline.
% 
% INPUTS
% data (statsdata):
% bondline_file (file dir): file containing 3 columns describing a 
%                           series of points on the bondplane.
%
%% -----------------------------------------------------------------------

% Get the data from the bondline file
bondline = importdata(bondline_file);

% Generate a best fit plane from the points and plot it to check.
[normal, ~, point] = affine_fit(bondline);


%% Plot Bestfit Plane
h = figure(1); clf(h); hold on; daspect([1 1 1]);

% plot the points
scatter3(bondline(:,1),bondline(:,2),bondline(:,3));

% plot a marker on plane and normal vector
plot3(point(1), point(2), point(3), 'ro', 'markersize', 5,...
      'markerfacecolor', 'red');
quiver3(point(1), point(2), point(3), 100*normal(1), 100*normal(2),...
        100*normal(3), 'r', 'linewidth', 2, 'AutoScale', 'off');

% plot the fitted plane
xmin = min(bondline(:,1));
xmax = max(bondline(:,1));
zmin = min(bondline(:,3));
zmax = max(bondline(:,3));

[X,Z] = meshgrid(linspace(xmin,xmax,3),linspace(zmin,zmax,3));

surf(X, -(normal(1)/normal(2)*X + normal(3)/normal(2)*Z -...
     dot(normal,point)/normal(2)), Z,'facecolor','red','facealpha',0.5);
%surf(X,repmat(point(2),size(X)),Z);

title('Best Fit Plane for Bondline'); hold off;

%% Calculate distances

[xyz,~,predu] = data.getPredictor({'x','y','z'});

[t,r,l] = reorient(xyz,normal,point);

data = data.addPredictor('t',t,predu{1});
data = data.addPredictor('b_dist',r,predu{1});
data = data.addPredictor('l',l,predu{1});
end

function [t,r,l] = reorient(xyz,normal,point)
% describe xyz in trl coordinates where the bondplane is the lt plane and
% r is the plane normal direction
fprintf('Calculating distances...');

% determine new basis interms of old basis
r_hat = normal/norm(normal); % already determined
l_hat = cross([1; 0; 0], r_hat);
l_hat = l_hat/norm(l_hat);
t_hat = cross(r_hat, l_hat);
t_hat = t_hat/norm(t_hat);

% compute change of basis rotation
R_10 = [t_hat, r_hat, l_hat];
R_01 = transpose(R_10);


% shift origin to bondline
xyz_shift = xyz - point;

% apply rotation to all points
trl =  R_01 * xyz_shift';
trl = trl';

t = trl(:,1) + min(trl(:,1));
r = trl(:,2);
l = trl(:,3) + min(trl(:,3));

fprintf(' DONE\n');
end