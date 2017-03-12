function [h] = pretty_bond_info(bond, key)
% Create a pretty bond information graphic showing adhesive penetration and
% bondline deformation together
%
% The figure consists of 4 panels:
% (A) L average of adhesive concentration in the RT plane
% (B) T average of adhesive concentration in the RL plane
% (C) L average of shear modulus in the RT plane
% (D) T average of shear modulus in the RL plane
%
%% Test input
% adhesive_map = random('unif', 0, 1, [200,300,500]) < 0.2;
% strain_map = random('unif', 0, 2e4, [200,300,500]);


%% Some parameters
bond = bond.get_bins();

ModLim = log10([1e-12, 1e-8]); % MPa

obuffer = 0.8; % in
ibuffer = 0.15; % in
side_lengths = (bond.limits(2,:) - bond.limits(1,:))...
               .* bond.resolution{1} * 2;

T = 1; R = 2; L = 3;


%%
h = figure();
A = subplot(2, 2, 1);
C = subplot(2, 2, 2);
B = subplot(2, 2, 3);
D = subplot(2, 2, 4);


%% A
subplot(A);
adhesive_fraction = permute(max(bond.adhesive_map, [],L, 'omitnan'),[T R L]);
adhesive_fraction(isnan(adhesive_fraction)) = 0;
contourf(bond.bins{R}, bond.bins{T}, adhesive_fraction);
%A.XDir = 'reverse';
A.XAxisLocation = 'Top';
colorbar('Location', 'northoutside');
caxis([0,1]);
% A.DataAspectRatio = [1, size(adhesive_map, T), 1];
title('Maximum Adhesive Volume Fraction', 'FontSize', 12);
ylabel(A, 'Tangential [mm]');
grid on;
grid minor;


%% B
subplot(B);
adhesive_fraction = permute(max(bond.adhesive_map, [],T, 'omitnan'),[L R T]);
adhesive_fraction(isnan(adhesive_fraction)) = 0;
contourf(bond.bins{R}, bond.bins{L}, adhesive_fraction);
%B.XDir = 'reverse';
B.YDir = 'reverse';
% colorbar('Location', 'westoutside');
caxis([0,1]);
% B.DataAspectRatio = [1, size(adhesive_map, T), 1];
xlabel(B, 'Radial [mm]');
ylabel(B, 'Longitudinal [mm]');
grid on;
grid minor;


%% C
subplot(C);
modulus_average = permute(mean(1./bond.modulus_map(key),L, 'omitnan'),[T R L]);
contourf(bond.bins{R}, bond.bins{T}, log10(abs(modulus_average * 1e-6)));
C.YAxisLocation = 'Right';
C.XAxisLocation = 'Top';
C.YDir = 'normal';
colorbar('Location', 'northoutside');
C.CLim = ModLim;
title(['Log Mean Compliance ', key ' Pa^{-1}'], 'FontSize', 12);
ylabel(C, 'Tangential [mm]');
grid on;
grid minor;


%% D
subplot(D);
modulus_average = permute(mean(1./bond.modulus_map(key),T, 'omitnan'),[L R T]);
contourf(bond.bins{R}, bond.bins{L}, log10(abs(modulus_average * 1e-6)));
D.YAxisLocation = 'Right';
D.YDir = 'reverse';
% colorbar('Location', 'eastoutside');
D.CLim = ModLim;
xlabel(D, 'Radial [mm]');
ylabel(D, 'Longitudinal [mm]');
grid on;
grid minor;


%%
B.Units = 'inches';
B.Position = [obuffer, obuffer, side_lengths(R), side_lengths(L)];

A.Units = 'inches';
A.Position = [obuffer, B.Position(2) + B.Position(4) + ibuffer,...
              side_lengths(R), side_lengths(T)];
          
D.Units = 'inches';
D.Position = [B.Position(1) + B.Position(3) + ibuffer, obuffer,...
              side_lengths(R), side_lengths(L)];

C.Units = 'inches';
C.Position = [B.Position(1) + B.Position(3) + ibuffer,...
              B.Position(2) + B.Position(4) + ibuffer,...
              side_lengths(R), side_lengths(T)];

h.Units = 'inches';
h.Position = [0, 0, 2*obuffer + ibuffer + 2*side_lengths(R),...
              1.5 + obuffer + ibuffer + side_lengths(T) + side_lengths(L)];

x = ibuffer/2 + side_lengths(R);
y = 1.25 + side_lengths(T);
text(A, x, y, bond.sample_name, 'Units', 'inches',...
     'HorizontalAlignment', 'center', 'FontSize', 14);
          
fprintf('PRETTY BOND: %s, %s\n', bond.sample_name, key);
end
