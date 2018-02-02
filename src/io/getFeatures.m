function data = getFeatures(info, data, window_size, scan_index)
% Calculate volume percentage of each feature in correlation window.
%
% PARAMETERS
% info - a json containing information about the sample
% data - the statsdata object
% window_size - the correlation window size
% scan_index - the index of this image in the series
% RETURNS
% data - altered statsdata object
%% -----------------------------------------------------------------------
left = data.vcount;
% For each subheading of "feature masks"
for feature = fieldnames(info.featureMasks)'
    feature_name = feature{1};

    try
        isdeformed = info.featureMasks.(feature_name).delta;
        k = scan_index;
    catch
        isdeformed = false;
        k = 1;
    end

    % Generate feature volumes from tagged slices.
    % Some tagged slices are in two parts (compound) and others are in one
    % (simple). 
    featureLabels = fieldnames(info.featureMasks.(feature_name));
    
    if ~strcmp(featureLabels{1}, 'image')
        
        fprintf(1, 'COMPOUND FEATURE: %s\n', feature_name);
        mask = logical([]);
        
        for i = 2:numel(featureLabels)
            
            mask = cat(3, mask, expand_compound(info, feature_name,...
                                                featureLabels{i}, k));
        end
        
    else
        
        fprintf(1, 'SIMPLE FEATURE: %s\n', feature_name);
        mask = expand_simple(info, feature_name, k);
    
    end

    % Reorder the indices if needed; TRL is the default order.
    unit = info.featureMasks.(feature_name).unit';
    mask = permute(mask,[...
                         strfind(cell2mat(unit), 'T'),...
                         strfind(cell2mat(unit), 'R'),...
                         strfind(cell2mat(unit), 'L'),...
                        ]);
                    
    assert(all(size(mask) == info.dimensions.value'),...
           'Mask size does not match volume dimensions');

    if isdeformed
        mask0 = expand_simple(info, feature_name, 1);
        mask0 = permute(mask0,[...
                               strfind(cell2mat(unit), 'T'),...
                               strfind(cell2mat(unit), 'R'),...
                               strfind(cell2mat(unit), 'L'),...
                              ]);
        data = calc_is_feature(data,feature_name,mask0,window_size,false);
        data = calc_is_feature(data,[feature_name num2str(k,'%02i')],...
                               mask,window_size, isdeformed);

        data.vcount = data.vcount - 1;
        initial = data.variables{data.vcount,2};
        final = data.variables{data.vcount+1,2};
        data.variables{data.vcount,2} = (final - initial) ./ initial;
        data.variables{data.vcount, 3} = 'percent';
        data.variables{data.vcount+1, 1} = [];
        data.variables{data.vcount+1, 2} = [];
        data.variables{data.vcount+1, 3} = [];

    else
        % place the things in the stuff.
        data = calc_is_feature(data,feature_name,mask,window_size, isdeformed);
    end
end
right = data.vcount;
[X, headers] = data.getPredictor(left+1:right);
figure, boxplot(X, headers);
title('Newly added features');
end

function mask = expand_simple(info, feature_name, scan_index)

% load and convert mask to boolean
mask = imstackload([info.dataLocation info.featureMasks.(feature_name).image]) > 0;

if ndims(mask) > 3
    
    channels = info.featureMasks.(feature_name).channel;
    
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
    
    mask = combined_channel;
end

mask = permute(mask,[2,1,3]);  % because MATLAB is different than ImageJ

if size(mask,3) > 1
    return;
else
    % slide it across the volume
    slide = info.featureMasks.(feature_name).slide';
    range = info.featureMasks.(feature_name).range';
    mask = shift(mask,slide,range);
end
end

function mask = expand_compound(info, feature_name, subfeature_name, scan_index)

% load and convert mask to boolean
mask = imstackload([info.dataLocation info.featureMasks.(feature_name).(subfeature_name).image]) > 0;

if ndims(mask) > 3
    
    channels = info.featureMasks.(feature_name).(subfeature_name).channel;
    
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
    
    mask = combined_channel;
end

mask = permute(mask,[2,1,3]);  % because MATLAB is different than ImageJ

if size(mask,3) > 1
    return;
else
    % slide it across the volume
    slide = info.featureMasks.(feature_name).(subfeature_name).slide';
    range = info.featureMasks.(feature_name).(subfeature_name).range';
    mask = shift(mask,slide,range);
end
end
