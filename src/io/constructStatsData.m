function constructStatsData(JSON, MAT_dir)
% Constructs a StatsData2 object for each of the time steps and each of the
% exports in the JSON.
%% -----------------------------------------------------------------------
disp(['CONSTRUCTING: ' JSON]);
info = readJSON(JSON);

for j = 2:numel(info.projectNumbers)
    scan = info.projectNumbers(j);

for i = 1:numel(info.digitalVolumeCorrelation.exports)
    export = info.digitalVolumeCorrelation.exports(i);

    
    %% LOAD THE DISPLACEMENT DATA FROM THE FILE
    datadir = sprintf('%s/recon_proj_%i/export.%04i', ...
                        info.dataLocation,...
                        scan,...
                        export);

    try
        data = loadVicVolume(datadir);
    catch
        fprintf('NOT FOUND: %s\n', datadir);
        continue
    end
    %% ADD CONSTANT LABELS
    try data = data.addLabel('time', info.experiment.time.value(j),...
                             info.experiment.time.unit);
    catch
    end
    try data = data.addLabel('humidity',...
                             info.experiment.humidity.value(j),...
                             info.experiment.humidity.unit);
    catch
    end
    try data = data.addLabel('cookTime',...
                             info.experiment.cookTime.value,...
                             info.experiment.cookTime.unit);
    catch
    end
    try data = data.addLabel('windowSize',...
                             info.digitalVolumeCorrelation.windowSize(i),...
                             'p');
    catch
    end
    try data = data.addLabel('force',...
                             info.experiment.force.value(j),...
                             info.experiment.force.unit);
    catch
    end
    try data = data.addLabel('species',...
                             info.species,...
                             '');
    catch
    end
    try data = data.addLabel('EP',...
                             info.bondline.EP,...
                             'p');
    catch
    end
    try data = data.addLabel('WP',...
                             info.bondline.WP,...
                             '1/p');
    catch
    end
    try data = data.addLabel('area',...
                             info.crossSection.area,...
                             info.crossSection.unit);
    catch
    end
    try data = data.addLabel('radial',...
                             info.crossSection.radial,...
                             info.crossSection.unit);
    catch
    end
    try data = data.addLabel('tangential',...
                             info.crossSection.tangential,...
                             info.crossSection.unit);
    catch
    end

    %% LABEL THE POINTS WITH FEATURES

    data = getFeatures(info, data, info.digitalVolumeCorrelation.windowSize(i), j);

    %% BONDLINE DISTANCE (PIXELS)
    try
        data = calc_bondline_distance(data, [info.dataLocation...
                                             info.bondline.points]);
        data = calc_penetration_half(data, info);
    catch
    end
    
    %% CALCULATE 3D MEAN STRAINS
    data = calc_mean_strains(data);
    
    %% CONVERT UNITS TO SI
    data = data.scale('x',info.resolution.value, info.resolution.unit);
    data = data.scale('y',info.resolution.value, info.resolution.unit);
    data = data.scale('z',info.resolution.value, info.resolution.unit);

    try
    data = data.scale('b_dist',info.resolution.value, info.resolution.unit);
    data = data.scale('t',info.resolution.value, info.resolution.unit);
    data = data.scale('l',info.resolution.value, info.resolution.unit);
    catch
    end
    
    try
    data = data.scale('time', 60, 's');
    catch
    end

    %% BOUNDING BOX DISTANCE (SI)
    try
        box = info.boundingBox .* info.resolution.value; % Txw,Ryh
        data = calc_distance_edge(data, box(1), box(2), box(3), box(4));
    catch
    end

    %% SAVE MODEL PARAMETERS
    datadir = sprintf('%s/recon_proj_%02i-%02i.mat', ...
                        MAT_dir,...
                        scan,...
                        export);
    if ~exist(MAT_dir,'dir')
        mkdir(MAT_dir)
    end
    save(datadir, 'data');
    fprintf(1, '%s\n', datadir);

    close('all')
end
end
