function h5d = plotCorrectSpiralDensity(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
%%
spiral_folder = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
frames_pre = 62:68; frames_post = 76:82;
hist_bin = 40;
framesN = numel(frames_pre);
load(fullfile(spiral_folder,'spiral_density_map_trial_types.mat'));
%%
labels_all = {'Correct pre-stim','Correct post-stim','Incorrect pre-stim','Incorrect post-stim',...
    'Miss pre-stim','Miss post-stim'};

h5d = figure('Renderer', 'painters', 'Position', [100 100 800 800]);
lineColor = 'k'; lineColor1 = 'w';
scale3 = 5/1;
hemi = [];
count = 1;
for j =1:3
    for k = 1:2
        if k==1
            unique_spirals = unique_spirals_pre{j};
        elseif k==2
            unique_spirals = unique_spirals_post{j};
        end
        ax1 = subplot(3,2,(j-1)*2+k);
        unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
        unique_spirals_unit = unique_spirals_unit./trialN(j)*35/framesN;                  % spirals/(mm^2*s)
        cmax(count,1) = max(unique_spirals_unit(:));
        scatter(unique_spirals(:,1),unique_spirals(:,2),...
            3,unique_spirals_unit,'filled');
        hold on;
        scale2=1;
        overlayOutlines(coords,scale2);
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
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
        title(labels_all{(j-1)*2+k});
        %%
        count = count+1;
    end
end

for j =1:6
    cmax2 = max(cmax);
    % cmax2 = cmax(1);
    ax1 = subplot(3,2,j);
    caxis([0,cmax2]);
    cb1 = colorbar;
    cb1.Ticks = [0,cmax2];
    cb1.TickLabels = {'0',num2str(round(cmax2*10)/10)};
end
%%
print(h5d, fullfile(save_folder,'Fig5d_task_spiral_density'), ...
    '-dpdf', '-bestfit', '-painters');