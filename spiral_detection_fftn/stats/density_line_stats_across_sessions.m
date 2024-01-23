% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
control_folder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\control';
permute_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\permute';
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
xx = 1:1140; yy = 1:1320;   
[xxq,yyq] = meshgrid(xx,yy);
count_sample_control = zeros(15,10,416);
count_sample_permute = zeros(15,10,416);
for kk = 1:session_total
    clear spiralsT
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    control = load(fullfile(control_folder,[mn '_density.mat']));
    permute = load(fullfile(permute_folder,[mn '_density.mat']));
    count = 1;
    for radius = 10:10:100
        clear unique_spirals_unit_control unique_spirals_unit_permute F_control F_permute vq1_control vq1_permute
        unique_spirals_unit_control = control.spiral_density{count}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = permute.spiral_density{count}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        % interp histgram counts 
        F_control = scatteredInterpolant(control.spiral_density{count}(:,1),control.spiral_density{count}(:,2),unique_spirals_unit_control);
        vq1_control = F_control(xxq,yyq);
        F_permute = scatteredInterpolant(permute.spiral_density{count}(:,1),permute.spiral_density{count}(:,2),unique_spirals_unit_permute);
        vq1_permute = F_permute(xxq,yyq);
        for i = 1:size(points,1)
            count_sample_control(kk,count,i) = vq1_control(points(i,2),points(i,1));
            count_sample_permute(kk,count,i) = vq1_permute(points(i,2),points(i,1));
        end
        count = count+1;
    end
end
%%
count_sample_control(count_sample_control<0) = 0;
count_sample_permute(count_sample_permute<0) = 0;
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
h = figure;
for radius = 1:10
    subplot(10,2,2*radius-1);
    for kk = 1:15
        plot(points_line,squeeze(count_sample_control(kk,radius,:)),'k');
        hold on;
        ylim([0,1]);
    end
    subplot(10,2,2*radius);
    for kk = 1:15
        plot(points_line,squeeze(count_sample_permute(kk,radius,:)),'k');
        hold on;
        ylim([0,1]);
    end
end
%%
max_control = max(count_sample_control,[],3);
max_permute = max(count_sample_permute,[],3);