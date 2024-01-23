%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%% 
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%%
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
%%
scale3 = 5;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
control_folder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\fftn_control\analysis';
permute_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\fftn_permutation\analysis';
load(fullfile(control_folder,'total_frames.mat'))
h = figure('Renderer', 'painters', 'Position', [100 100 900 900]);
count = 1;
for radius = 10:10:50
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
    ax1 = subplot(4,5,count);
    scatter(control.unique_spirals(:,1),control.unique_spirals(:,2),3,unique_spirals_unit_control,'filled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax1,hot);
    axis image; axis off;
    xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
    caxis([0,cmax]);
    cb = colorbar;    
    cb.Ticks = [0,cmax];
    cb.TickLabels = {'0',num2str(round(cmax*10)/10)};
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)+0.05 a(2) a(3) a(4)/2])% To change size
    %%
    ax2 = subplot(4,5,5+count);
    scatter(permute.unique_spirals(:,1),permute.unique_spirals(:,2),3,unique_spirals_unit_permute,'filled');
    hold on;
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax2,hot);
    axis image; axis off;
    xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
    caxis([0,cmax])
    cb = colorbar; 
    cb.Ticks = [0,cmax];
    cb.TickLabels = {'0',num2str(round(cmax*10)/10)};
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)+0.05 a(2) a(3) a(4)/2])% To change size
    count = count+1;
end
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
    ax1 = subplot(4,5,10+count);
    scatter(control.unique_spirals(:,1),control.unique_spirals(:,2),3,unique_spirals_unit_control,'filled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax1,hot);
    axis image; axis off;
    xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
    caxis([0,cmax]);
%     xtick = 0:200:1000; xtick_string = string(xtick);
%     xticks(xtick); xticklabels(xtick_string);
%     ytick = 0:200:1200; ytick_string = string(ytick);
%     yticks(ytick); yticklabels(ytick_string);
    cb = colorbar; 
    cb.Ticks = [0,cmax];
    cb.TickLabels = {'0',num2str(round(cmax*10)/10)};
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)+0.05 a(2) a(3) a(4)/2])% To change size
    %%
    ax2 = subplot(4,5,15+count);
    scatter(permute.unique_spirals(:,1),permute.unique_spirals(:,2),3,unique_spirals_unit_permute,'filled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax2,hot);
    axis image; axis off;
    xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
    caxis([0,cmax])
%     xtick = 0:200:1000; xtick_string = string(xtick);
%     xticks(xtick); xticklabels(xtick_string);
%     ytick = 0:200:1200; ytick_string = string(ytick);
%     yticks(ytick); yticklabels(ytick_string);
    cb = colorbar; 
    cb.Ticks = [0,cmax];
    cb.TickLabels = {'0',num2str(round(cmax*10)/10)};
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)+0.05 a(2) a(3) a(4)/2])% To change size
    count = count+1;
end
%%
print(h, 'spiral_fftn_radius_all2', '-dpdf', '-bestfit', '-painters');