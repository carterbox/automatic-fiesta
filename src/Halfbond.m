classdef Halfbond < Bond
% An object for storing half bond information
%
%
%% -----------------------------------------------------------------------
    properties
        % Specimen information
        is_top = NaN;
        
        modulus_map = [];
        adhesive_map = [];
        cell_density_map = [];
    end
    
    methods
        function self = Halfbond(bond, is_top)
             
            self.sample_name = bond.sample_name;
            self.species = bond.species;
            self.resolution = bond.resolution;
            self.units = bond.units;
            self.is_top = is_top;
            
            % Tranformation from global to TRL coordinates
            self.T = bond.T;
            self.shift = bond.shift;

            % Data in TRL coordinates
            self.strains = bond.strains;
            if self.is_top
                self.adhesive = bond.adhesive(bond.adhesive <= 0);
                self.cellwall = bond.cellwall(bond.cellwall <= 0);
                
                for i = 1:numel(self.strains)
                    try
                    keep = self.strains{i}.getData('r') <= 0;
                    self.strains{i} = self.strains{i}.reduce(keep);
                    catch
                        fprintf(1,'Failed to reduce strain %i.\n', i);
                    end
                end
            else
                self.adhesive = bond.adhesive(bond.adhesive >= 0);
                self.cellwall = bond.cellwall(bond.cellwall >= 0);
                
                for i = 1:numel(self.strains)
                    try
                    keep = self.strains{i}.getData('r') >= 0;
                    self.strains{i} = self.strains{i}.reduce(keep);
                    catch
                        fprintf(1,'Failed to reduce strain %i.\n', i);
                    end
                end
            end
        end
        
        
        function save(self, datadir)
            
            if self.is_top
                tag = '_lo';
            else
                tag = '_hi';
            end
            
            halfbond = self;
            datadir = [datadir '/' self.sample_name tag];
            
            save(datadir, 'halfbond', '-v7.3');
        end
        
        
        function map = get_adhesive_map(self)
            if ~isempty(self.adhesive_map)
                map = self.adhesive_map;
                return
            end
            
            
        end
        
        
        function map = get_cell_density_map(self)
            if isempty(self.cell_density_map)
                map = self.cell_density_map;
                return
            end
        end
        
        
        function map = get_modulus_map(self)
            if isempty(self.modulus_map)
                map = self.modulus_map;
                return
            end
        end
    end

end

function [G] = bulk_shear_modulus(data, Tyz0)

m1 = fitlm([data.getData('Eyz'); zeros(data.ocount, 1)],...
           [data.getData('Tyz'); ones(data.ocount, 1) * Tyz0],...
           'VarNames', {'Eyz', 'Tyz'});
G = m1.Coefficients.Estimate(2);
% m1.plotSlice

G = abs(G);
end

function [Gmap] = local_shear_modulus(data, Tyz0)
grid_size = [0.05, 0.05, 0.05]; % mm


points = data.getData({'x', 'b_dist', 'z'});
points(:,2) = abs(points(:,2));

[gridx, gridy, gridz] = ndgrid(min(points(:,1)):grid_size(1):max(points(:,1)),...
                               min(points(:,2)):grid_size(2):max(points(:,2)),...
                               min(points(:,3)):grid_size(3):max(points(:,3)));

Gmap = zeros(size(gridx));
parfor i = 1:numel(gridx)
    
    local = data.reduce(gridx(i) <= points(:,1) & points(:,1) < gridx(i)+grid_size(1) ...
                      & gridy(i) <= points(:,2) & points(:,2) < gridy(i)+grid_size(2) ...
                      & gridz(i) <= points(:,3) & points(:,3) < gridz(i)+grid_size(3));
    if local.ocount > 0
        Gmap(i) = bulk_shear_modulus(local, Tyz0);
    end
end
end