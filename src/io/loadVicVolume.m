function [ sdata ] = loadVicVolume( directory )
%LOADVICVOLUME returns a StatsData2 object of the data in an export
% directory from vic volume.

%% -----------------------------------------------------------------------
if ~exist(directory), error('%s does not exist!', directory); end

% Loads the names of all the files in the directory
fcontents = dir(directory);
addpath( genpath(directory) );
start = 1;
while fcontents(start).isdir == 1
    start = start + 1;
    if start > length(fcontents), error('There are no files there!'); end
end

fprintf('READ: %s\n', directory);
count = 0; data = [];
for i = start:length(fcontents)
    %checks that the entry is not a directory
    if(~fcontents(i).isdir)
        fprintf('%s ',fcontents(i).name);
        if rem(i-2,5) == 0; fprintf('\n'); end
        count = count + 1;
        temp = importdata(fcontents(i).name);
        data = cat(1, data, temp.data);
    end
end
fprintf('LOADED %i FILES\n', count);

%% Separate bad correlations. Status will be non-zero
headers = {'σ [voxel]','ε_xx','ε_yy','ε_zz','ε_xy','ε_xz','ε_yz'};
figure(); boxplot(data(:,7:13), 'Labels',headers);
title('Raw VicVolume data');
data = data(data(:,7)~=-1,:);

% for i = 7:6 %7:13
%     q = prctile(data(:,i),[5 95]);
%     bad = data(:,i) < q(1) | data(:,i) > q(2);
%     data = data(~bad,:);
%     fprintf('REMOVE OUTLIARS: %s, %0.2f%%\n', headers{i-6},sum(bad)/numel(bad)*100);
% end
figure(); boxplot(data(:,7:13), 'Labels',headers);
title('Reduced VicVolume data');

%% First filling in missing points

% if ~isempty(missing_points)
% parfor f = 8:size(data,2)x y z u v w sigma Exx Eyy
% Ezz Exy Exz Eyz.
%     d0 = data(:,f);
%     
%     % Interpolate the values of F at the missing points
%     F = griddatan(x,dg(:,f-7),missing_points);
% 
%     j = 1;
%     for missing_row = find(bad)'
%         d0(missing_row) = F(j);
%         j = j + 1;
%     end
%     
%     data(:,f) = d0;
% end
% end
% 
% fprintf(' FILLED %i POINTS\n',numel(bad));

sdata = statsdata2();
sdata = sdata.addPredictor('x',data(:,1),'p');
sdata = sdata.addPredictor('y',data(:,2),'p');
sdata = sdata.addPredictor('z',data(:,3),'p');

sdata = sdata.addResponse('u',data(:,4),'p');
sdata = sdata.addResponse('v',data(:,5),'p');
sdata = sdata.addResponse('w',data(:,6),'p');

sdata = sdata.addResponse('sigma',data(:,7),'');
sdata = sdata.addResponse('Exx',data(:,8),'');
sdata = sdata.addResponse('Eyy',data(:,9),'');
sdata = sdata.addResponse('Ezz',data(:,10),'');
sdata = sdata.addResponse('Exy',data(:,11),'');
sdata = sdata.addResponse('Exz',data(:,12),'');
sdata = sdata.addResponse('Eyz',data(:,13),'');
end
