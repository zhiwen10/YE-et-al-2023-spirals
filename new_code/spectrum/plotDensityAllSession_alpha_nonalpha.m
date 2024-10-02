%%
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_fft_bump');
filteredSpirals_alpha_all = [];
filteredSpirals_nonalpha_all = [];
alpha_frames = 0;
nonalpha_frames = 0;
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    subfolder = [mn '_' tdb '_' num2str(en)];
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,[fname '_histogram.mat']));
    filteredSpirals_alpha_all = [filteredSpirals_alpha_all; filteredSpirals_alpha];
    filteredSpirals_nonalpha_all = [filteredSpirals_nonalpha_all; filteredSpirals_nonalpha];
    alpha_frames = alpha_frames+numel(alphaFrames);
    nonalpha_frames = nonalpha_frames+numel(nonalphaFrames);
end
%%
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_fft_bump_density\allSessions');
hist_bin = 40;
label = 'alpha';
frameN = alpha_frames;
for radius = 40:10:100
    clear unique_spirals
    spirals_temp = [];
    spirals_temp = filteredSpirals_alpha_all(filteredSpirals_alpha_all(:,3)==radius,:);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
    fname = ['histogram_radius_' num2str(radius) '.mat'];
    save(fullfile(save_folder,label,fname),...
        'unique_spirals','frameN');
end

label = 'nonalpha';
frameN = nonalpha_frames;
for radius = 40:10:100
    clear unique_spirals
    spirals_temp = [];
    spirals_temp = filteredSpirals_nonalpha_all(filteredSpirals_nonalpha_all(:,3)==radius,:);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
    fname = ['histogram_radius_' num2str(radius) '.mat'];
    save(fullfile(save_folder,label,fname),...
        'unique_spirals','frameN');
end
%%
local_data_folder = fullfile(data_folder,...
    'spirals\spirals_freq\spirals_fft_bump_density\allSessions');
alpha_data_folder = fullfile(local_data_folder,'alpha');
nonalpha_data_folder = fullfile(local_data_folder,'nonalpha');
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
%% get cortical atlas path and tree
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
xx = 1:1140; yy = 1:1320;
[xxq,yyq] = meshgrid(xx,yy);
[maskPath,st] = get_cortex_atlas_path(data_folder);
root1 = '/997/';
ctx = '/997/8/567/688/';
scale3 = 5;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
hs3c = figure('Renderer', 'painters', 'Position', [100 100 1100 600]);
count = 1;
for radius = 40:10:100
    alpha = load(fullfile(alpha_data_folder,['histogram_radius_' num2str(radius) '.mat']));
    nonalpha = load(fullfile(nonalpha_data_folder,['histogram_radius_' num2str(radius) '.mat']));
    %%
    unique_spirals_unit_alpha = alpha.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_alpha = unique_spirals_unit_alpha./frameN*35; % spirals/(mm^2*s)\
    %% interp histgram counts 
    Fna = scatteredInterpolant(alpha.unique_spirals(:,1),...
        alpha.unique_spirals(:,2),...
        unique_spirals_unit_alpha);
    vq1_na = Fna(xxq,yyq);
    for i = 1:size(points,1)
        count_sample_alpha_all(i,count) = vq1_na(points(i,2),points(i,1));
    end    
    %%   
    unique_spirals_unit_nonalpha = nonalpha.unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit_nonalpha = unique_spirals_unit_nonalpha./frameN*35; % spirals/(mm^2*s)
    %% interp histgram counts 
    Fna = scatteredInterpolant(nonalpha.unique_spirals(:,1),...
        nonalpha.unique_spirals(:,2),...
        unique_spirals_unit_nonalpha);
    vq1_na = Fna(xxq,yyq);
    for i = 1:size(points,1)
        count_sample_nonalpha_all(i,count) = vq1_na(points(i,2),points(i,1));
    end  
    %%
    cmax = max([unique_spirals_unit_alpha;unique_spirals_unit_nonalpha]);
    %%
    ax1 = subplot(3,7,count);
    scatter(alpha.unique_spirals(:,1),alpha.unique_spirals(:,2),3,unique_spirals_unit_alpha,'filled');
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
    set(cb,'Position',[a(1)+0.035 a(2) a(3) a(4)/2])% To change size
    %%
    ax2 = subplot(3,7,7+count);
    scatter(nonalpha.unique_spirals(:,1),nonalpha.unique_spirals(:,2),3,unique_spirals_unit_nonalpha,'filled');
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
    set(cb,'Position',[a(1)+0.035 a(2) a(3) a(4)/2])% To change size
    %%
    subplot(3,7,14+count);
    plot(points_line,count_sample_alpha_all(:,count),'r');
    hold on;
    plot(points_line,count_sample_nonalpha_all(:,count),'k');
    xlim([0,800]);
    ylim([0,0.3]);
    yticks([0,0.1,0.2,0.3]);
    yticklabels({'0' '0.1' '0.2' '0.3'});
    xticks([0 200 400 600 800]);
    xticklabels({'0','2','4','6','8'});
    ylabel('sprials/(mm^2*s)');  
    %%
    count = count+1;
end
print(hs3c, fullfile(save_folder,'desnity_map_alpha_nonalpha'),...
    '-dpdf', '-bestfit', '-painters');