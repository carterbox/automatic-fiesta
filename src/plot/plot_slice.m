function plot_slice(info, MAT_dir, modeltag, predictors)

% generate a background image for the plot
stack_dir = sprintf('%s/recon_proj_%i', info.dataLocation, info.projectNumbers(1)); 
back_stack = imstackload(stack_dir);
back_gray = uint16(mean(back_stack,3));
back_rad = uint16(mean(back_stack,2));
back_color = repmat(back_gray, [1,1,3]);
rad_color = repmat(squeeze(back_rad)', [1,1,3]);
clear back_stack back_gray

%% ITERATE THROUGH ALL THE DATA
k = 1; DATA = statsdata();
for j = 2:numel(info.projectNumbers)
    scan = info.projectNumbers(j);

for i = 1:numel(info.digitalVolumeCorrelation.exports)
    export = info.digitalVolumeCorrelation.exports(i);
    
    datadir = sprintf('%s/recon_proj_%02i-%02i.mat',...
                        MAT_dir,...
                        scan,...
                        export);
    disp(k);
    load(datadir, 'data');
    k = k + 1;

%% This is a section
[T, ~, ~] = data.getPredictor('x');
[R, ~, ~] = data.getPredictor('y');
[TQ,RQ] = ndgrid(min(T)*1000/1.1:5:max(T)*1000/1.1,...
                min(R)*1000/1.1:5:max(R)*1000/1.1);

fraction = numel(TQ) / data.ocount / 30;
RADIUS = numel(TQ)/ 3000;


strains = {'Em', 'Exx', 'Eyy'};
for k = 1:numel(strains)
RESPONSE = strains{k};
TITLE = sprintf('%s--(N%i-W%i-T%03i)--%s',...
                info.sample,...
                export,...
                info.digitalVolumeCorrelation.windowSize(i),...
                info.experiment.time.value(j),...
                RESPONSE);

% LOAD A MODEL FOR COMPARISON
load(['./models/',RESPONSE,'_',modeltag,'.mat']);
P = model.predict(data.getPredictor(predictors));
% =========

[S, ~, ~] = data.getResponse(RESPONSE);
x = random('bino', 1, fraction, size(S)) > 0;

drawnow;

% 1 =========
fig_boxplot = figure(1); clf; % plot a boxplot to show the range of strains.
if exist('P','var')
    boxplot([S,P],{'DVC','Model'});
else
    boxplot(S);
end
title(TITLE);
% ylim([-0.1,0.1]);

drawnow;
saveas(fig_boxplot,[info.dataLocation '/' TITLE '_boxplot'],'png');

% 4 =========
% plot a 2D averaged contour of ALL data over the average image.
Si = scatteredInterpolant(T(x)*1000/1.1,R(x)*1000/1.1,S(x)); 
Pi = scatteredInterpolant(T(x)*1000/1.1,R(x)*1000/1.1,P(x)); 

SQ = Si(TQ,RQ); 
PQ = Pi(TQ,RQ); 

fig_2DContour = figure(4); clf;
subplot(2,1,1);
imshow(back_color); hold on;

cmax = 0.15;
colormap parula
colorbar();
caxis([-cmax,0]);
[c,h] = contourf(TQ,RQ,SQ);
ch = get(h,'child'); alpha(ch,0.1)


title(TITLE);
drawnow;
subplot(2,1,2);
imshow(back_color); hold on;

cmax = 0.15;
colormap parula
colorbar();
caxis([-cmax,0]);
[c,h] = contourf(TQ,RQ,PQ);
ch = get(h,'child'); alpha(ch,0.1)


title(modeltag);
drawnow;
% tightfig;
saveas(fig_2DContour,[info.dataLocation '/' TITLE '_2DContour'],'png');

% % 2 =========
% % plot a 2D scatter of random subset on top of the average image. 
% fig_random2DScatter = figure(2); clf;
% imshow(back_color); hold on;
% 
% colormap parula
% colorbar();
% caxis([-cmax,cmax]);
% scatter(T(x)*1000/1.1,R(x)*1000/1.1,RADIUS, S(x),'filled');
% 
% title(TITLE);
% 
% hold off;
% 
% drawnow;
% saveas(fig_random2DScatter,[info.dataLocation '/' TITLE '_random2DScatter'],'png');

% 3 =========
% plot a 3D scatter of random subset with the average image on the z plane
fig_random3DScatter = figure(3); clf;
% subplot(2,1,1);

scatter3(T(x)*1000/1.1,R(x)*1000/1.1,S(x),RADIUS,S(x),'filled'); hold on;
scatter3(T(x)*1000/1.1,R(x)*1000/1.1,P(x),RADIUS,[0 0 0],'filled'); hold off;

zlim([-size(rad_color,1)/10000,0]);
ylim([0,size(rad_color,2)]);
view([1 0 0]);
xlabel('tangential');
ylabel('radial');
title(TITLE);
pbaspect([1000 size(rad_color,2) size(rad_color,1)])
colormap parula
% colorbar();
caxis([-cmax,0]);

% tightfig
drawnow;

% subplot(2,1,2);
% imshow(rad_color);
% drawnow;

saveas(fig_random3DScatter,[info.dataLocation '/' TITLE '_random3DScatter'],'png');
savefig(fig_random3DScatter,[info.dataLocation '/' TITLE '_random3DScatter']);

disp('END OF LOOP');
end
end
end

close all;

end
