function hs15j = plotCorrectSpiralPrePost(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
% frames_pre = 59:68; frames_post = 75:84;
% frames_pre = 64:66; frames_post = 78:80;
frames_pre = 63:67; frames_post = 77:81;
framesN = numel(frames_pre);
labels = {'Correct pre-stim','Correct post-stim'};
hist_bin = 40;
%%
load(fullfile(data_folder,'task','spirals','CorrectSpiralsPrePost.mat'));
%%
hs15j = figure('Renderer', 'painters', 'Position', [100 100 800 800]);
for j =1:2
    unique_spirals = unique_spirals_all{j};
    ax1 = subplot(1,2,j);
    unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./trialN(j)*35/framesN;                 % spirals/(mm^2*s)
    cmax(j,1) = max(unique_spirals_unit(:));
    scatter(unique_spirals(:,1),unique_spirals(:,2),...
        3,unique_spirals_unit,'filled');
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse');
    colormap(ax1,hot);
    % axis off; 
    axis image;
    xlim(ax1,[0,1140]);
    ylim(ax1,[0,1320]);
    xtick = 0:200:1000;
    xtick_string = string(xtick);
    xticks(xtick);
    xticklabels(xtick_string);
    ytick = 0:200:1200;
    ytick_string = string(ytick);
    yticks(ytick);
    yticklabels(ytick_string);
    title(labels{j});
end

for j =1:2
    cmax2 = max(cmax);
    ax1 = subplot(1,2,j);
    caxis([0,cmax2]);
    cb1 = colorbar;
    cb1.Ticks = [0,cmax2];
    cb1.TickLabels = {'0',num2str(round(cmax2*10)/10)};
end
%%
print(hs15j, fullfile(save_folder,'FigS15j_task_spirals_pre_post'), ...
    '-dpdf', '-bestfit', '-painters');