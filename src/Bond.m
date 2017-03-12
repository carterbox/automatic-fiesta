classdef Bond
% An object for storing bond information
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
        grid_spacing = [0.05, 0.05, 0.05]; % mm
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
        
        modulus_map = containers.Map;
        strain_map = containers.Map;
        adhesive_map = [];
        cell_density_map = [];
        
        % Load transfer data
        model = NaN;
    end
    
    methods
        function self = Bond(JSON)
            self.JSON = JSON;
            self = self.add_info();
            self.modulus_map = containers.Map;
            self.strain_map = containers.Map;
        end
        
        function self = add_info(self)
            info = readJSON(self.JSON);
            
            self.sample_name = info.sample;
            self.species = info.species;
            self.resolution = {info.resolution.value,...
                               info.resolution.unit};
            
            try
                self.force = info.experiment.force.value;
                self.tangential = {info.crossSection.tangential,...
                                   info.crossSection.unit};
            catch
            end
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
        
        function self = add_adhesive(self)
        % Return trl coordinates of adhesive labeled pixels
            fprintf(1, 'ADD ADHESIVE\n');
            xyz = get_feature_coordinates(self.JSON, 'adhesive');
            
            [trl, self] = self.xyz_to_trl(xyz);
            
            trl = self.within_limits(trl);
            
            self.adhesive = trl;
        end
        
        function self = add_cellwall(self)
            % Return the trl coordinates of cell wall labeled pixels
            
            xyz = get_feature_coordinates(self.JSON, 'cellwall');
            
            [trl, self] = self.xyz_to_trl(xyz);

            trl = self.within_limits(trl);
            
            self.cellwall = trl;
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
                    if isa(scan, 'cell')
                        datadir = sprintf('%s/%s/export.%04i', ...
                                            info.dataLocation,...
                                            scan{1},...
                                            export);
                        scan = scan{1};
                    else
                        datadir = sprintf('%s/recon_proj_%i/export.%04i', ...
                                            info.dataLocation,...
                                            scan,...
                                            export);
                        scan = num2str(scan);
                    end

                    try
                        data = loadVicVolume(datadir);
                    catch
                        fprintf('NOT FOUND: %s\n', datadir);
                        continue
                    end
                    
                    data = self.rotate_strains(data);
                    
                    data = calc_von_mises(data);
                    
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
                    fprintf(1, 'NEW RANGES: %s | %s\n', self.sample_name, scan);
                    fprintf(1,'TRL MIN: %.2f, %.2f, %.2f\n', new_min);
                    fprintf(1,'TRL MAX: %.2f, %.2f, %.2f\n', new_max);
                    
                    self.limits(1,:) = min(new_min, self.limits(1,:));
                    self.limits(2,:) = max(new_max, self.limits(2,:));
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
        
        function self = get_strain_map(self, key)
            % Bin strain data into grid for a given strain
            
            if ~self.strain_map.isKey(key)
                
                fprintf(1, 'BIN STRAIN: %s\n', key);
                
                strain_map = [];

                grid_spacing = self.grid_spacing;
                
                if strcmp('u', key) || strcmp('v', key) || strcmp('w', key)
                 grid_spacing = grid_spacing * 4;
                end
                
                for i = 2:numel(self.strains)
                    
                    % collect DVC data for the requested 
                    trl = self.strains{i}.getData({'t','r','l'});
                    E_ = self.strains{i}.getData(key);
                    
                    [~, ~, ~, map] = binned_data(trl.*self.resolution{1},...
                                                 E_,...
                                                 self.limits.*self.resolution{1},...
                                                 grid_spacing, @mean);
                                             
                    % save the result
                    strain_map = cat(4, strain_map, map);
                end
                
                self.strain_map(key) = strain_map;
            end
            
        end
        
        function self = get_modulus_map(self, key)
            % Calculate modulus map from strain map
            
            if ~self.strain_map.isKey(key)
                self = self.get_strain_map(key);
            end
            
            if ~self.modulus_map.isKey(key)
                fprintf(1, 'MODULUS MAP: %s\n', key);
            
                % compute applied engineering shear stress
                Tyz = self.force ./ (self.tangential{1} * 5.0); % N/mm^2 == MPa
                
                nansize = size(self.strain_map(key));
                modulus_map = nan(nansize(1:3));
                sorted_strain = permute(self.strain_map(key), [4,1,2,3]);
                
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
                
                self.modulus_map(key) = modulus_map;
            end
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
        
        function self = savefig(self, dir, key)
            % saves a pretty bond info
            
            self = self.get_modulus_map(key);
            
            h = pretty_bond_info(self, key);
            saveas(h, [dir '/' self.sample_name], 'png');
        end
        
        
        function trl = within_limits(self, trl)
            % reduce trl data to only include points within the limits
            
            keep = prod(trl >= self.limits(1,:), 2)...
                 & prod(trl <= self.limits(2,:), 2);
             
            trl = trl(keep, :);
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
        
        function [xyz] = trl_to_xyz(self, trl)
            xyz = (self.T' * (trl - [self.shift(1), 0, self.shift(3)])')' +...
                   self.shift;
        end
        
        
        function [data, self] = rotate_strains(self, data)
            
            % calculate bondline if not already
            if isempty(self.T) || isempty(self.shift)
                self = self.add_bondline();
            end

            Exx = data.getResponse('Exx');
            Eyy = data.getResponse('Eyy');
            Ezz = data.getResponse('Ezz');
            Exy = data.getResponse('Exy');
            Exz = data.getResponse('Exz');
            Eyz = data.getResponse('Eyz');

            parfor i = 1:numel(Exx)
    
                S = [Exx(i),Exy(i),Exz(i);Exy(i),Eyy(i),Eyz(i);Exz(i),Eyz(i),Ezz(i)]; 

                Sp = self.T * S;
                
                Exx(i) = Sp(1,1);
                Exy(i) = Sp(1,2);
                Exz(i) = Sp(1,3);
                Eyy(i) = Sp(2,2);
                Eyz(i) = Sp(2,3);
                Ezz(i) = Sp(3,3);
            end

            data = data.addResponse('Exx',Exx,'');
            data = data.addResponse('Eyy',Eyy,'');
            data = data.addResponse('Ezz',Ezz,'');
            data = data.addResponse('Exy',Exy,'');
            data = data.addResponse('Exz',Exz,'');
            data = data.addResponse('Eyz',Eyz,'');


            [uvw, ~, units] = data.getResponse({'u', 'v', 'w'});
            
            uvw1 = (self.T * uvw')';
            
            data = data.addResponse('u', uvw1(:, 1), units{1});
            data = data.addResponse('v', uvw1(:, 2), units{2});
            data = data.addResponse('w', uvw1(:, 3), units{3});

        end
        
        
        function get_slices(self, location)
            
            % calculate the middle slice of the region
            lmean = mean(self.limits, 1);
            
            % explicitly list the coordinates of the points in these slices
            [trt, trr, trl] = ndgrid(self.limits(1,1):1:self.limits(2,1),...
                        self.limits(1,2):1:self.limits(2,2),...
                        lmean(3));
            
            tr = cat(2, trt(:), trr(:), trl(:));
            
            [rlt, rlr, rll] = ndgrid(lmean(1),...
                        self.limits(1,2):1:self.limits(2,2),...
                        self.limits(1,3):1:self.limits(2,3));
                    
            rl = cat(2, rlt(:), rlr(:), rll(:));
           
                    
            % convert these coordinates to xyz
            
            xy = int16(self.trl_to_xyz(tr));
            
            yz = int16(self.trl_to_xyz(rl));
            
            % get the coordinates from the volumes
            
            volume = imstackload(location);
            volume = permute(volume, [2,1,3]);
            
            xy(xy < 1) = 1;
            xy(xy(:,1) > size(volume,1), 1) = size(volume, 1);
            xy(xy(:,2) > size(volume,2), 2) = size(volume, 2);
            xy(xy(:,3) > size(volume,3), 3) = size(volume, 3);
            % save the slices as images
            i = sub2ind(size(volume), xy(:,1), xy(:,2), xy(:,3));
            
            imgxy = permute(reshape(volume(i), size(trt)), [1,2]);
            
            yz(yz < 1) = 1;
            yz(yz(:,1) > size(volume,1), 1) = size(volume, 1);
            yz(yz(:,2) > size(volume,2), 2) = size(volume, 2);
            yz(yz(:,3) > size(volume,3), 3) = size(volume, 3);
            % save the slices as images
            j = sub2ind(size(volume), yz(:,1), yz(:,2), yz(:,3));
            
            imgyz = permute(reshape(volume(j), size(rlt)), [3,2,1]);
            
            imwrite(imgxy,[self.sample_name '_tr.png']);
            imwrite(imgyz,[self.sample_name '_rl.png']);
        end
        
        function [EP, WP, lEP, lWP] = penetration_stats(self)
            % Return EP and WP metrics

            area = prod(self.limits(2,[1,3]) - self.limits(1,[1,3]));
            distance = abs(self.adhesive(:, 2)); % radial component
            t_coord = self.adhesive(:, 1); % distance from edge
            t_space = self.grid_spacing(1) / self.resolution{1};
            lEP = zeros(size(min(t_coord):t_space:max(t_coord)));
            lWP = zeros(size(min(t_coord):t_space:max(t_coord)));
            
            for i = 1:size(self.bins{1}) 
                lo = min(t_coord) + (i-1) * t_space;
                hi = lo + t_space;
                x = t_coord > lo & t_coord < hi;
                larea = (self.limits(2,3) - self.limits(1,3)) * t_space;
                lEP(i) = effective_penetration(distance(x), larea) * self.resolution{1};
                lWP(i) = weighted_penetration(distance(x)) * self.resolution{1};
                i = i + 1;
            end
            
%             histogram(lEP); hold on;
%             histogram(lWP);
            EP = mean(lEP);
%             EPu = std(lEP);
            WP = mean(lWP);
%             WPu = std(lWP);
            lWP = squeeze(lWP);
            lEP = squeeze(lEP);
            
%             EP = effective_penetration(distance, area) * self.resolution{1};
%             WP = weighted_penetration(distance) * self.resolution{1};
        end

        function [ds, Eyzs] = plot_shear_drop(self, key, N, centers, Lrange, dir, step)
           % plots the shear strain as a function of distance from the
           % notch
           
           h = figure(); hold on;
           mark = {'d-','s-','o-'};
           self = self.get_bins();
           ds = {};
           Eyzs = {};
           leg = {};
           
            
           for i = 1:numel(centers)
               % Choose the two bins closest to the bondline
               Rbins = find(self.bins{2} > centers(i), 1) - 1;
               Rbins = Rbins - N:Rbins + N;

               % Average across T direction
               Eyz = self.strain_map(key);
               Eyz = Eyz(:, Rbins, :, step);
               Eyz = squeeze(mean(mean(Eyz, 1, 'omitnan'), 2, 'omitnan'));

               d = self.bins{3};
               keep = d >= Lrange(1) & d <= Lrange(2);
               d = d(keep);
               Eyz = 1./Eyz(keep) * 1e-6;

               plot(d, Eyz, mark{self.species});
               ds{i} = d;
               Eyzs{i} = Eyz;
               leg{i} = num2str(centers(i));
           end
           title([self.sample_name ' load transfer']);
           xlabel('notch distance [mm]');
           ylabel(['mean ' key]);
           legend(leg)
           hold off;
           
           saveas(h, [dir '/' self.sample_name '.drop.png'], 'png');
           
           ds = cell2mat(ds);
           Eyzs = cell2mat(Eyzs);
        end
    
    end
end

%% Utilites

function EP = effective_penetration(distance, area)
% Calculate the volume of pixels per bondline area [m^3/m^2]

fprintf('\nCalculating EP...');

EP = numel(distance)/area;

fprintf(' %f', EP);
end

function WP = weighted_penetration(distance)
% Calculate the weighted penetration of the points in the volume

fprintf('\nCalculating WP...');

WP = sqrt( sum(distance.^2)/ numel(distance) ); % [sqrt(m^2 m^3/m^3)]

fprintf(' %f\n', WP);
end

function G = estimate_modulus(stress, strain)
% Return the estimated slope of a line given by stress strain

if numel(strain) > 1
    m1 = fitlm(strain, stress, 'VarNames', {'Epsilon', 'Sigma'});
    G = m1.Coefficients.Estimate(2);
%     m1.plotSlice
%     G = abs(G);
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

