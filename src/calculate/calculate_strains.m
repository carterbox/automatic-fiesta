function calculate_strains(d, strain_file)
%CALCULATE_STRAINS takes the matrix d which contains locations and
% displacement data for a regular grid.
%
% version 0.1.1
% d: a matrix whose column headings are:
% 1, 2, 3, 4, 5     , 6     , 7, 8, 9,...
% n, x, y, z, status, objmin, u, v, w,...(other values)
%
% strain_file: the name of the output file to put all the data in.
%
%% Determine the shape of the grid
d = d(:,1:9);
d = sortrows(d,[4 3 2]);

a = (unique(d(:,2)));
b = (unique(d(:,3)));
c = (unique(d(:,4)));

h1 = mean(diff(a));
h2 = mean(diff(b));
h3 = mean(diff(c));

a = numel(a);
b = numel(b);
c = numel(c);

parfor i = 1:3
    UVW{i} = reshape(d(:,i+6),a,b,c);
end

%% Then smooth the data
fprintf('Smooth data...');

parfor i = 1:3
    UVW{i} = smooth3(UVW{i},'gaussian');
end

fprintf(' DONE\n');
%% Differentiate the Regular Data
fprintf('Differentiate...');

parfor i = 1:3
    [Dx{i},Dy{i},Dz{i}] = gradient(UVW{i},h1,h2,h3);
end

clear UVW
fprintf(' DONE\n');
%% Convert data to Tensor

Vx = Dx{2}; Wx = Dx{3}; Ux = Dx{1};
Vy = Dy{2}; Wy = Dy{3}; Uy = Dy{1};
Vz = Dz{2}; Wz = Dz{3}; Uz = Dz{1};
clear Dx Dy Dz

S = cat(2,Ux(:),Vy(:),Wz(:),(Vz(:)+Wy(:))./2,(Uz(:)+Wx(:))./2,(Uy(:)+Vx(:))./2);

% figure(1), quiver3(d(:,1),d(:,2),d(:,3),Ux(:),Vy(:),Wz(:)); daspect([1,1,1]);
% figure(2), histogram(Ux(:)); axis([-1/2 1/2 -inf inf]);
% figure(3), histogram(Vy(:)); axis([-1/2 1/2 -inf inf]);
% figure(4), histogram(Wz(:)); axis([-1/2 1/2 -inf inf]);

clear Ux Uy Uz Vx Vy Vz Wx Wy Wz

%% Save data to file

results = cat(2, d(:,1:4), S);
results = sortrows(results,1);
results = results';

file = fopen(strain_file,'w','n','UTF-8');

fprintf(file, '% 6s % 8s % 8s % 8s % 10s % 10s % 10s % 10s % 10s % 10s\r\n','n','x','y','z','Exx','Eyy','Ezz','Eyz','Exz','Exy');
fprintf(file, '% 6i % 8.2f % 8.2f % 8.2f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f\r\n',results);

fclose(file);

end
