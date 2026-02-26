function hs3c = plotMapDataVsFftn(data_folder,save_folder,freq)
%%
freq_name = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
%%
local_data_folder = fullfile(data_folder,...
    'spirals','spirals_fftn');
control_data_folder = fullfile(local_data_folder,...
    freq_name,'control_map');
fftn_data_folder = fullfile(local_data_folder,...
    freq_name,'fftn_map');
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
%% get cortical atlas path and tree
[maskPath,st] = get_cortex_atlas_path(data_folder);
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
scale3 = 5;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
load(fullfile(control_data_folder,'total_frames.mat'))
hs3c = figure('Renderer', 'painters', 'Position', [100 100 900 900]);
count = 1;
for radius = 10:10:100
    control = load(fullfile(control_data_folder,['histogram_radius_' num2str(radius) '.mat']));
    permute = load(fullfile(fftn_data_folder,['histogram_radius_' num2str(radius) '.mat']));
    %%
    unique_spirals_unit_control = control.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_control = unique_spirals_unit_control./frame_all*35; % spirals/(mm^2*s)
    unique_spirals_unit_permute = permute.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_permute = unique_spirals_unit_permute./frame_all*35; % spirals/(mm^2*s)
    %%
    cmax = max([unique_spirals_unit_control;unique_spirals_unit_permute]);
    %%
    if radius>= 60
        count1 = count+5;
    else
        count1 = count;
    end
    ax1 = subplot(4,5,count1);
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
    ax2 = subplot(4,5,5+count1);
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
%%
save_name = fullfile(save_folder,...
    ['FigS3c_spirals_radius_' freq_name '.pdf']);
print(hs3c, save_name,'-dpdf', '-bestfit', '-painters');
end