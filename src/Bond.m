classdef Bond
% An object for storing bond information
%
%
%% -----------------------------------------------------------------------
    properties
        % Specimen information
        species = NaN;
        sample_name = '';
        file_location = '';
        JSON = '';
        
        % spatial information
        resolution = NaN; % mm
        grid_spacing = [0.050, 0.050, 0.050]; % mm
        tangential = NaN; % mm
        force = []; % N
        
        % Tranformation from global to TRL coordinates
        T = [];
        shift = [];
        
        % Data in TRL coordinates
        limits = [];
        strains = {};
        adhesive = [];
        cellwall = [];
        bins = {};
        
        modulus_map = [];
        strain_map = [];
        adhesive_map = [];
        cell_density_map = [];
    end
    
    methods
        function self = Bond(JSON)
            self.JSON = JSON;
            self = self.add_info();
        end
        
        function self = add_info(self)
            info = readJSON(self.JSON);
            
            self.sample_name = info.sample;
            self.species = info.species;
            self.resolution = {info.resolution.value,...
                               info.resolution.unit};
            
            self.force = info.experiment.force.value;
            self.tangential = {info.crossSection.tangential,...
                               info.crossSection.unit};
        end
        
        function self = add_strains(self)
            % Load strain data and set trl limits
            fprintf(1, 'ADD STRAINS\n');
            
            info = readJSON(self.JSON);
            
            self.limits = nan(2, 3);
            
            % iterate through each project and export
            for j = 2:numel(info.projectNumbers)
                scan = info.projectNumbers(j);
                for i = 1:numel(info.digitalVolumeCorrelation.exports)
                    export = info.digitalVolumeCorrelation.exports(i);

                    % LOAD THE DISPLACEMENT DATA FROM THE FILE
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
                    
                    [xyz, ~, units] = data.getData({'x','y','z'});
                        
                    [trl, self]= self.xyz_to_trl(xyz);
                    
%                     figure, boxplot(xyz);
%                     figure, boxplot(trl);

                    data = data.addData('t', trl(:, 1), units{1});
                    data = data.addData('r', trl(:, 2), units{1});
                    data = data.addData('l', trl(:, 3), units{1});

                    self.strains{j} = data;
                    
                    % adjust the limits to include this data
                    new_min = min(trl,[],1);
                    new_max = max(trl,[],1);
                    fprintf(1, 'NEW RANGES: %s | %i\n', self.sample_name, scan);
                    fprintf(1,'TRL MIN: %.2f, %.2f, %.2f\n', new_min);
                    fprintf(1,'TRL MAX: %.2f, %.2f, %.2f\n', new_max);
                    
                    self.limits(1,:) = min(new_min, self.limits(1,:));
                    self.limits(2,:) = max(new_max, self.limits(2,:));
                end
            end
        end
        
        function self = get_strain_map(self, name)
            % Bin strain data into grid for a given strain
            fprintf(1, 'BIN STRAIN: %s\n', name);
            
            if isempty(self.strain_map)
                for i = 2:numel(self.strains)
                    
                    % collect DVC data for the requested 
                    trl = self.strains{i}.getData({'t','r','l'});
                    E_ = self.strains{i}.getData(name);
                    
                    [~, ~, ~, map] = binned_data(trl.*self.resolution{1},...
                                                 E_,...
                                                 self.limits.*self.resolution{1},...
                                                 self.grid_spacing, @mean);
                                             
                    % save the result
                    self.strain_map = cat(4, self.strain_map, map);
                end
            end
            
        end
        
        function self = get_bins(self)
           % calculate the centers of the bins
           
           limits = self.limits.*self.resolution{1};
           grid_spacing = self.grid_spacing;
           
           % add half grid spacing to put locations at centers of bins
           self.bins{1} = symmetric_range(limits(1,1),grid_spacing(1),limits(2,1)) + grid_spacing(1)/2;
           self.bins{2} = symmetric_range(limits(1,2),grid_spacing(2),limits(2,2)) + grid_spacing(2)/2;
           self.bins{3} = symmetric_range(limits(1,3),grid_spacing(3),limits(2,3)) + grid_spacing(3)/2;
        end
        
        function self = get_modulus_map(self)
            % Calculate modulus map from strain map
            fprintf(1, 'MODULUS MAP\n');
            
            if isempty(self.strain_map)
                self = self.get_strain_map('Eyz');
            end
            
            if isempty(self.modulus_map)
            
                % compute applied engineering shear stress
                Tyz = self.force ./ (self.tangential{1} * 5.0); % N/mm^2 == MPa
                
                modulus_map = nan(size(self.strain_map(:,:,:,1)));
                sorted_strain = permute(self.strain_map, [4,1,2,3]);
                
                parfor i = 1:numel(modulus_map)
                    strain = sorted_strain(:, i);
                    
                    % preppend no strain for reference step
                    strain = [0; strain];
                    
                    % remove steps/locations that have no strain value
                    keep = ~isnan(strain);
                    strain = strain(keep);
                    stress = Tyz(keep);
                    
                    modulus_map(i) = estimate_modulus(stress, strain);
                end
                
                self.modulus_map = modulus_map;
            end
        end
        
        function self = add_adhesive(self)
        % Return trl coordinates of adhesive labeled pixels
            fprintf(1, 'ADD ADHESIVE\n');
            xyz = get_feature_coordinates(self.JSON, 'adhesive');
            
            [trl, self] = self.xyz_to_trl(xyz);
            
            trl = self.within_limits(trl);
            
            self.adhesive = trl;
        end
        
        function self = get_adhesive_map(self)
            % Bin adhesive data into grid
            fprintf(1, 'BIN ADHESIVE\n');
            
            if isempty(self.adhesive_map)
                
                [~, ~, ~, map] = binned_data(self.adhesive.*self.resolution{1},...
                                             1,...
                                             self.limits.*self.resolution{1},...
                                             self.grid_spacing, @sum);
                
                % convert from pixel counts to volume density              
                map = map .* self.resolution{1}^3 ./ prod(self.grid_spacing); 
                
                % save the result
                self.adhesive_map = map;
            end
        end
        
        function self = add_cellwall(self)
            % Return the trl coordinates of cell wall labeled pixels
            
            xyz = get_feature_coordinates(self.JSON, 'cellwall');
            
            [trl, self] = self.xyz_to_trl(xyz);

            trl = self.within_limits(trl);
            
            self.cellwall = trl;
        end
        
        function trl = within_limits(self, trl)
            % reduce trl data to only include points within the limits
            
            keep = prod(trl >= self.limits(1,:), 2)...
                 & prod(trl <= self.limits(2,:), 2);
             
            trl = trl(keep, :);
        end
        
        function self = add_bondline(self)
            % Compute transformation for reorienting coordinates
            
            info = readJSON(self.JSON);
            bondline_file = [info.dataLocation info.bondline.points];
            
            % Get the data from the bondline file
            bondline = importdata(bondline_file);

            % Generate a best fit plane from the points and plot it to check.
            [normal, ~, point] = affine_fit(bondline);

            plot_best_plane(bondline, point, normal);

            % Prevent directions from flipping
            if normal(2) < 0
                normal = -normal;
            end
            
            % Determine new basis interms of old basis
            r_hat = normal/norm(normal); % already determined
            l_hat = cross([1; 0; 0], r_hat);
            l_hat = l_hat/norm(l_hat);
            t_hat = cross(r_hat, l_hat);
            t_hat = t_hat/norm(t_hat);

            % Compute change of basis rotation
            R_10 = [t_hat, r_hat, l_hat];
            R_01 = transpose(R_10);
        
            self.T = R_01;
            self.shift = point;
        end
        
        function save(self, datadir)
            % Save the Bond to a file.
            
            if isempty(self.file_location)
                datadir = [datadir '/' self.sample_name];
                self.file_location = datadir;
            else
                datadir = self.file_location;
            end
            
            bond = self;
            save(datadir, 'bond', '-v7.3');
        end
        
        function [trl, self] = xyz_to_trl(self, xyz)
           % convert xyz coordinates to trl coordinates
            
           % calculate bondline if not already
            if isempty(self.T) || isempty(self.shift)
                self = self.add_bondline();
            end

           trl = (self.T * (xyz - self.shift)')'...
                        + [self.shift(1), 0, self.shift(3)];
        end
        
        function savefig(self, dir)
            % saves a pretty bond info
            h = pretty_bond_info(self);
            saveas(h, [dir '/' self.sample_name], 'png');
        end
    end

    

end

function G = estimate_modulus(stress, strain)
% Return the estimated slope of a line given by stress strain

if numel(strain) > 1
    m1 = fitlm(strain, stress, 'VarNames', {'Epsilon', 'Sigma'});
    G = m1.Coefficients.Estimate(2);
    % m1.plotSlice
    G = abs(G);
else
    G = NaN;
end

end

function [xbins, ybins, zbins, A] = binned_data(xyz, value, limits, grid_spacing, func)
% Put an Nx3 dataset into a binned grid based on the limits 2x3 and the
% grid spacing 3x1. func determines how the values in the bins are
% combined. value is the value of the data at each point.

% There's a better way to make the bins which involves doing it after
% accumulation.
xbins = symmetric_range(limits(1,1),grid_spacing(1),limits(2,1));
ybins = symmetric_range(limits(1,2),grid_spacing(2),limits(2,2));
zbins = symmetric_range(limits(1,3),grid_spacing(3),limits(2,3));
                           
binned_xyz = floor(xyz ./ grid_spacing);
binned_xyz = binned_xyz - floor(limits(1, :) ./ grid_spacing) + 1;

A = accumarray(binned_xyz, value,...
               [numel(xbins), numel(ybins), numel(zbins)], func, NaN);

end

function r = symmetric_range(lo, step, hi)
% Return a range from lo to high where zero is the center.

r = (floor(lo/step)*step:step:floor(hi/step)*step)';

end

function xyz = get_feature_coordinates(JSON, feature_name)
    info = readJSON(JSON);
    
    fprintf(1, 'SIMPLE FEATURE: %s\n', feature_name);
    mask = expand_simple(info, feature_name, 1);

    % Reorder the indices if needed; TRL is the default order.
    unit = info.featureMasks.(feature_name).unit';
    mask = permute(mask,[...
                         strfind(cell2mat(unit), 'T'),...
                         strfind(cell2mat(unit), 'R'),...
                         strfind(cell2mat(unit), 'L'),...
                        ]);
                    
    assert(all(size(mask) == info.dimensions.value'),...
           'Mask size does not match volume dimensions');
    
    [x,y,z] = ind2sub(size(mask), find(mask));
    
    xyz = [x,y,z];
end

function mask = expand_simple(info, feature_name, scan_index)
% load and convert mask to boolean

mask = imstackload([info.dataLocation info.featureMasks.(feature_name).image{scan_index}]) > 0;

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

function plot_best_plane(bondline, point, normal)
    h = figure(1); clf(h); hold on; daspect([1 1 1]);
    
    % plot the points
    a0 = scatter3(bondline(:,1),bondline(:,2),bondline(:,3));

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
    
    a0.Parent.YDir = 'reverse';
    title('Best Fit Plane for Bondline'); hold off;
end

