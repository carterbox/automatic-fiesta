function data = calc_penetration_half(data, info)
%CALC_PENETRATION adhesive penetration characterizer.
%
%INPUTS
% VOLUME (boolean): A boolean volume where true values are the locations of
% adhesive.
% POINTS (Nx3 array): Points in the volume that approximate the bondline.
%
%OUTPUTS
% EP (double) Effective Penetration: A surface density; the average mass
% per unit bondline.
% WP (double) Weighted Penetration: Similar to the second moment of area;
% accounts for both mass and perpendicular distance. Masses farther away
% from the bondline count more.
%
% version 0.0.1 - fixed calculations so they operate on XZ planes, and
% fixed graphics so they show at the proper scale.
%% TEST CASE

% points = [1,1,1;0,0,1;1,0,1];
% volume = true(3,3,3);

% volume = imstackload('/media/OCT14B/OCT14B/Subsets/Daniel/HPPHEL00/Link to color_01');
% volume = volume(:,:,:,1) | volume(:,:,:,3);
% points = importdata('/media/OCT14B/OCT14B/Subsets/Daniel/HPPHEL00/bondline.csv');

%% ------------------------------------------------------------------------

mask = imstackload([info.dataLocation info.featureMasks.adhesive.image{1}]);

channels = info.featureMasks.adhesive.channel;
    
combined_channel = zeros(size(mask,1),...
                         size(mask,2), size(mask,3), 'logical' );

if contains(channels, 'R')
    fprintf('ADD CHANNEL RED\n')
    combined_channel = combined_channel | mask(:,:,:,1);
end
if contains(channels, 'G')
    fprintf('ADD CHANNEL GREEN\n')
    combined_channel = combined_channel | mask(:,:,:,2);
end
if contains(channels, 'B')
    fprintf('ADD CHANNEL BLUE\n')
    combined_channel = combined_channel | mask(:,:,:,3);
end

volume = combined_channel; clear mask combined_channel

points = importdata([info.dataLocation info.bondline.points]);

%% CALCULATE a best fit plane

[normal, ~, point] = affine_fit(points);

%% PLOT BEST FIT PLANE
h1 = figure(); hold on;
title('Best Fit Plane for Bondline');
daspect([1 1 1]);

% points
scatter3(points(:,1),points(:,2),points(:,3));

% normal vector
plot3(point(1), point(2), point(3), 'ro', 'markersize', 5, 'markerfacecolor', 'red');
quiver3(point(1), point(2), point(3), 100*normal(1), 100*normal(2), 100*normal(3), 'r','linewidth',2, 'AutoScale','off');

% fitted plane
xmin = min(points(:,1));
xmax = max(points(:,1));
zmin = min(points(:,3));
zmax = max(points(:,3));

[X,Z] = meshgrid(linspace(xmin,xmax,3),linspace(zmin,zmax,3));

surf(X, -(normal(1)/normal(2)*X + normal(3)/normal(2)*Z -...
     dot(normal,point)/normal(2)), Z,'facecolor','red','facealpha',0.5);

drawnow;
hold off;

%% CALCULATE EP & WP

% For each TRUE value in volume calculate the WP and WP
[EP, WP, density] = penetration(volume, normal, point);


%%

d = data.getData('b_dist');

ep = double(d < 0) .* EP(1) + double(d >= 0) .* EP(2);
wp = double(d < 0) .* WP(1) + double(d >= 0) .* WP(2);

data = data.addData('EP', ep, 'p');
data = data.addData('WP', wp, 'p');
data.other{1} = density;

end

function [EP, WP, density] = penetration(volume, normal, point)
% Calculate the effective and weighted penetration of the points in the
% volume for each side of the plane defined by a normal and a point.
%
%% BONDPLANE AREA
% Calculate the area of the plane. The plane normal contains a z
% component. Otherwise, this calculated area will be infinite.
if(normal(2) == 0), error('Bondline should horizontal in slices!' + ...
                          ' Caclulated area is infinite!'); end

% Switch x,y because matlab imageJ coordinate difference
[~,xmax,zmax] = size(volume);
area = sqrt((-normal(1)/normal(2)).^2 + (-normal(3)/normal(2)).^2 + 1) * (xmax - 1) * (zmax - 1);

%% ADHESIVE DISTANCES
trus = find(volume);
% Switch x,y because matlab imageJ coordinate difference
[y,x,z] = ind2sub(size(volume),trus);
clear trus;

normal = normal/norm(normal);
distance = (x-point(1))*normal(1) + (y-point(2))*normal(2) + (z-point(3))*normal(3);

%% EFFECTIVE PENETRATION
% The average amount of adhesive per unit of bondplane. Can be added
% directly. sum( Vi ) / BondArea
EP = [0, 0];
EP(1) = sum(distance <  0)/area;
EP(2) = sum(distance >= 0)/area;

fprintf('CALCULATE EP: %f, %f\n', EP(1), EP(2)); % [m^3/m^2]

%% WEIGHTED PENETRATION
% The 2-root of the weighted average squared distance of adhesive from
% the bondline. Cannot be added directly. Because they are weighted
% averages. sqrt( sum( yi^2 * Vi ) / sum( Vi ) )
WP = [0, 0];
WP(1) = sqrt( sum(distance(distance <  0).^2)/ sum(distance <  0) );
WP(2) = sqrt( sum(distance(distance >= 0).^2)/ sum(distance >= 0) );

fprintf('CALCULATE WP: %f, %f\n', WP(1), WP(2)); % [m^2 * m^3/m^3]

%% DENSITY PLOT
density = zeros(xmax,zmax,2);

N = 10;
ysteps = floor((0:N)*(size(volume,1))/N);
for i = 1:numel(ysteps)-1
    clear x y z
    [y,x,z] = ndgrid(ysteps(i)+1:ysteps(i+1),1:size(volume,2),1:size(volume,3));
    d = (x-point(1))*normal(1) + (y-point(2))*normal(2) + (z-point(3))*normal(3);
    density(:,:,1) = density(:,:,1) + squeeze(sum(volume(ysteps(i)+1:ysteps(i+1),:,:) & d <  0 , 1));
    density(:,:,2) = density(:,:,2) + squeeze(sum(volume(ysteps(i)+1:ysteps(i+1),:,:) & d >= 0 , 1));
end

%% PLOT THINGS
h2 = figure();
histogram(distance);
title('Particle Distance from Bondline');
drawnow;

h3 = figure();
scatter(EP, WP);
daspect([1,1,1]);
title('Effective Penetration vs Weighted Penetration');
xlabel('EP [m^3/m^2]');
ylabel('WP [m^2/m^3]');

h4 = figure();

subplot(1,2,1);
title('Adhesive Density (-)');
contourf(permute(density(:,:,1),[2,1]));
set(gca,'xaxislocation','top','ydir','reverse');
daspect([1,1,1]);
xlabel('Tangential');
ylabel('Longitudinal');

subplot(1,2,2);
title('Adhesive Density (+)');
contourf(permute(density(:,:,2),[2,1]));
set(gca,'xaxislocation','top','ydir','reverse');
daspect([1,1,1]);
xlabel('Tangential');
ylabel('Longitudinal');

end