function [ J ] = shift( I, slide, range)
%Generates a stack of images, J, such that slice located at z moves along
%the line defined by (x0,y0,z0) and (x1,y1,z1) in the range z = [z_min, z_max).

cslide = num2cell(slide);
crange = num2cell(range);
[x0, y0, z0, x1, y1, z1, z_] = cslide{:};
[z_min, z_max] = crange{:};

dx = (x1-x0)/(z1-z0);
dy = (y1-y0)/(z1-z0);

J = zeros([size(I), z_max-z_min], 'like', I);

for z = z_min:(z_max-1)
    d = floor([dx,dy]*(z-z_));
    J(:,:,z-z_min+1) = circshift(I,d);
end

end

