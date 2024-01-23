%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
control_folder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\fftn_control\analysis';
permute_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\fftn_permutation\analysis';
load(fullfile(control_folder,'total_frames.mat'))
h = figure('Renderer', 'painters', 'Position', [100 100 600 900]);
count = 1;
for radius = 60:10:100
    control = load(fullfile(control_folder,['histogram_control_radius_' num2str(radius) '.mat']));
    permute = load(fullfile(permute_folder,['histogram_fftn_radius' num2str(radius) '.mat']));
    %%
    unique_spirals_unit_control = control.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_control = unique_spirals_unit_control./frame_all*35; % spirals/(mm^2*s)
    unique_spirals_unit_permute = permute.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_permute = unique_spirals_unit_permute./frame_all*35; % spirals/(mm^2*s)
    %%
    cmax = max([unique_spirals_unit_control;unique_spirals_unit_permute]);
    %%
    ax1 = subplot(5,2,2*count-1);
    scatter(control.unique_spirals(:,1),control.unique_spirals(:,2),3,unique_spirals_unit_control,'filled');
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    colormap(ax1,hot);
    axis image; axis off;
    xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
    caxis([0,cmax]);
    xtick = 0:200:1000; xtick_string = string(xtick);
    xticks(xtick); xticklabels(xtick_string);
    ytick = 0:200:1200; ytick_string = string(ytick);
    yticks(ytick); yticklabels(ytick_string);
    colorbar;    
    %%
    ax2 = subplot(5,2,2*count);
    scatter(permute.unique_spirals(:,1),permute.unique_spirals(:,2),3,unique_spirals_unit_permute,'filled');
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    colormap(ax2,hot);
    axis image; axis off;
    xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
    caxis([0,cmax])
    xtick = 0:200:1000; xtick_string = string(xtick);
    xticks(xtick); xticklabels(xtick_string);
    ytick = 0:200:1200; ytick_string = string(ytick);
    yticks(ytick); yticklabels(ytick_string);
    colorbar;
    count = count+1;
end
%%
print(h, 'spiral_fftn_radius_60_100', '-dpdf', '-bestfit', '-painters');