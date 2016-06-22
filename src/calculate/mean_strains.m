
% data = loadvicvolume('data_0001.csv');
% 
% Exx = data(:,8); Eyy = data(:,9); Ezz = data(:,10);
% Eyz = data(:,13); Exy = data(:,11); Exz = data(:,12);

function [mean_strain, Exx, Eyy, Ezz] = mean_strains(Exx,Eyy,Ezz,Eyz,Exy,Exz,num_axis)
%MEAN_STRAINS takes the 6 components of the 3x3 symmetric strain tensor and
% returns the mean strain and the 3 principle strains. The mean strain can
% be demoted to a biaxial mean strain by including the optional parameter 
% num_axis.
% version 0.0.0
%% -----------------------------------------------------------------------
if nargin < 7; num_axis = 3; end

mean_strain = zeros(size(Exx));
parfor i = 1:numel(Exx)
    e = eig([Exx(i),Exy(i),Exz(i);Exy(i),Eyy(i),Eyz(i);Exz(i),Eyz(i),Ezz(i)]); 
    mean_strain(i) = mean(e(1:num_axis));
    Exx(i) = e(1);
    Eyy(i) = e(2);
    Ezz(i) = e(3);
end

clear e i
end

