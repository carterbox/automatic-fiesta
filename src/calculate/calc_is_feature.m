function data = calc_is_feature(data, label, map, bsize, isdeformed)
% Determine the volume fraction of the correlation window that belongs to
% a feature.
% PARAMETERS
% data -    statsdata object which contains the x,y,z locations of the points
%           of interest.
% label -   the name of the label for the volume fractions when they are
%           added to the statsdata object
% map -     a feature map of the feature values are [0, 1]
% bsize -   the size of the region around the point to include in the volume
%           fraction calculation.
% isdeformed -  determines whether deformations are applied to the
%               volumes before they go to the map to grab the labels.

%% Labels points with booleans on whether they are located on a feature. Uses
% map to decide whether each point exists on the feature. Label is the name of
% the feature.

if nargin < 5
    isdeformed = false;
end

x = data.getPredictor('x');
y = data.getPredictor('y');
z = data.getPredictor('z');

if ~isdeformed && bsize > 1
    % region is square and has volume greater than 1
    
    % define region around each point
    lo = floor(bsize/2) - 1;
    hi = ceil(bsize/2);
    [sR, sC, sP, eR, eC, eP] = deal(x-lo,y-lo,z-lo,x+hi,y+hi,z+hi);
    
    regionSum = -ones(size(sR));
    parfor i = 1:numel(sR)
        region = map(sR(i):eR(i),sC(i):eC(i),sP(i):eP(i));
        regionSum(i) = sum(region(:));
    end
    assert(all(regionSum >= 0));
    
elseif false
    % region is non-square
    % Apply displacements and strains.
    
    % define relative region around each point
    lo = -floor(bsize/2); hi = -lo;
    [dx, dy, dz] = ndgrid(lo:hi,lo:hi,lo:hi);
    dx = dx(:); dy = dy(:); dz = dz(:);
    
    % get deformations for each region
    u = data.getResponse({'u', 'v', 'w'});

    Exx = data.getResponse('Exx');
    Eyy = data.getResponse('Eyy');
    Ezz = data.getResponse('Ezz');
    Exy = data.getResponse('Exy');
    Exz = data.getResponse('Exz');
    Eyz = data.getResponse('Eyz');
    
    regionSum = -ones(size(x));
    parfor i = 1:data.ocount
        
%         E = [Exx(i), Exy(i), Exz(i);...
%              Exy(i), Eyy(i), Eyz(i);...
%              Exz(i), Eyz(i), Ezz(i)];
        
        t = repmat(u(i,:) + [x(i), y(i), z(i)], numel(dx), 1);
        
        % deform each subregion
%         x1 = [dx, dy, dz]*E + t;
        x1 = [dx, dy, dz] + t;
        
        % round to nearest integer location
        x1 = round(x1);
        
        % use 1D indicies of 3D locations
        a = uint32(sub2ind(size(map),x1(:,1),x1(:,2),x1(:,3)));
        region = map(a);
        regionSum(i) = sum(region(:));
    end
    assert(all(regionSum >= 0));
    
elseif isdeformed
    % Apply displacements only.
    
    % get deformations for each region
    u = data.getResponse({'u', 'v', 'w'});
    u = round(u);
    
    % define region around each point
    lo = floor(bsize/2) - 1;
    hi = ceil(bsize/2);
    [sR, sC, sP, eR, eC, eP] = deal(x - lo + u(:,1),...
                                    y - lo + u(:,2),...
                                    z - lo + u(:,3),...
                                    x + hi + u(:,1),...
                                    y + hi + u(:,2),...
                                    z + hi + u(:,3));
    
    regionSum = -ones(size(sR));
    parfor i = 1:numel(sR)
        region = map(sR(i):eR(i),sC(i):eC(i),sP(i):eP(i));
        regionSum(i) = sum(region(:));
    end
    assert(all(regionSum >= 0));
    
else
    % each point is of size 1
    a = uint32(sub2ind(size(map),x,y,z));
    regionSum = map(a);
end

fraction = regionSum / (bsize^3);
assert(all(fraction <= 1));
    
fprintf(1,'%f percent of points are on %s.\n',...
        sum(fraction(:)) / numel(fraction) * 100, label);

data = data.addPredictor(label, fraction, 'bool');
end

function [J] = integralimage3D(V)
% INTEGRALIMAGE3D computes the integral image of a 3D volume. Which is an
% array whose indecies contain the sum of all the intensities in the
% original image that are above and to the left. This one is inclusive and
% 1 indexed.
% 
% INPUT
% V: the volume.
%
% OUTPUT
% J (double): the integral image of V.
%
% NOTES
% http://en.wikipedia.org/wiki/Summed_area_table
% We need the extra bit depth for large volumes.
% Example: 128*1024^3 ~ 1.3x10^11 > 2^32
% Also, Matlab processes double faster than int64.
%% -----------------------------------------------------------------------

% Pad the integral image with zeros to make the loop simpler (There's no
% need for special cases at the edges). Make the default value -1 for
% troubleshooting purposes.
[x0,y0,z0] = size(V);
%J = padarray(-ones(x0,y0,z0, 'double'), [1,1,1], 'pre');
J(x0+1,y0+1,z0+1) = double(0);
s = size(J);
J = tall(J);

V = V;

for i = 2:x0+1
for j = 2:y0+1
for k = 2:z0+1
    % Single pass integral image computation extrapolated to 3D from
    % Wikipedia.
    J(i,j,k) = V(i-1,j-1,k-1) + J(sub2ind(s,i-1,j,k)) + J(sub2ind(s,i,j-1,k)) + J(sub2ind(s,i,j,k-1))...
    - J(sub2ind(s,i,j-1,k-1)) - J(sub2ind(s,i-1,j,k-1)) - J(sub2ind(s,i-1,j-1,k)) + J(sub2ind(s,i-1,j-1,k-1));
end
end
end
end