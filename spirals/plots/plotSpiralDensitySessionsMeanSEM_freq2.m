function hs3h = plotSpiralDensitySessionsMeanSEM_freq2(freq,data_folder,save_folder)
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
load(fullfile(data_folder,'spirals','spirals_freq','spirals_density_line',...
    ['spiralDensityLine_' freq_folder,'.mat']))
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
load(fullfile(data_folder,'spirals','spirals_density',...
    'histogram_40pixels.mat'),'frame_total');                              % load total frame N calculated previously
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
mean_g7 = mean(max_density(1:11));
std_g7 = std(max_density(1:11));
mean_g8 = mean(max_density(12:15));
std_g8 = std(max_density(12:15));
[hh,p,ci,stats] = ttest2(max_density(1:11),max_density(12:15));
%%
session_total = 15;
mean_spirals = mean(count_sample,2);
std_spirals = std(count_sample,[],2)./sqrt(session_total);
%%
% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
%%
hs3h = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
shadedErrorBar(points_line(:,1),mean_spirals,std_spirals,'lineProps','g');
xlim([0,800]);
ylim([0,3]);
yticks([0,0.5,1,1.5,2,2.5,3]);
yticklabels({'0' '0.5' '1' '1.5' '2' '2.5' '3'});
xticks([0 200 400 600 800]);
xticklabels({'0','2','4','6','8'});
ylabel('sprials/(mm^2*s)');  
title(freq_folder,'interpreter','none');
print(hs3h, fullfile(save_folder,['FigS3h_spiralsDensityLine_' freq_folder '.pdf']), ...
    '-dpdf', '-bestfit', '-painters');