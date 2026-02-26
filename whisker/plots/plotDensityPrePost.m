function [h5d,h5f] = plotDensityPrePost(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
scale = 1;
[row,col] = find(BW);
brain_index = [col,row];
%%
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
load(fullfile(data_folder,'whisker','spirals_peri_stim','Whisker_spirals_pre_post.mat'));
%% right SSp index
clear areaPath
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
Utransformed = projectedAtlas1;
scale = 1;

hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp_right = [col,row];

hemi = 'left';
[indexSSp2,UselectedSSp2] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexSSp2);
indexSSp_left = [col2,row2];
%%
% areaPath(1) = "/997/8/567/688/695/315/500/985/"; %MOp
areaPath(1) = "/997/8/567/688/695/315/500/993/"; % MOs
MOsArea = strcat(areaPath(:));

% sensoryArea = [];
Utransformed = projectedAtlas1;
scale = 1;
hemi = 'right';
[indexMOs,UselectedSSp] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexMOs);
indexMOs_right = [col,row];

hemi = 'left';
[indexMOs,UselectedSSp2] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexMOs);
indexMOs_left = [col2,row2];
%%
areaPath(1) = "/997/8/567/688/695/315/500/985/"; %MOp
MOsArea = strcat(areaPath(:));
% sensoryArea = [];
Utransformed = projectedAtlas1;
scale = 1;
hemi = 'right';
[indexMOp,UselectedSSp] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexMOp);
indexMOp_right = [col,row];

hemi = 'left';
[indexMOp,UselectedSSp2] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexMOp);
indexMOp_left = [col2,row2];
%%
for kk = 1:5
    spirals_post_temp = spirals_post_cell{kk,1};
    [lia1,locb1] = ismember(spirals_post_temp(:,1:2),indexSSp_right,'rows');
    spirals_post_temp2 = spirals_post_temp(lia1,:);
    ratio_post(kk,1) = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
    %%
    spirals_pre_temp = spirals_pre_cell{kk,1};
    [lia1,locb1] = ismember(spirals_pre_temp(:,1:2),indexSSp_right,'rows');
    spirals_pre_temp2 = spirals_pre_temp(lia1,:);
    ratio_pre(kk,1) = sum(spirals_pre_temp2(:,4) == -1)./size(spirals_pre_temp2,1);
end
for kk = 1:5
    spirals_post_temp = spirals_post_cell{kk,1};
    [lia1,locb1] = ismember(spirals_post_temp(:,1:2),indexSSp_left,'rows');
    spirals_post_temp2 = spirals_post_temp(lia1,:);
    ratio_post2(kk,1) = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
    %%
    spirals_pre_temp = spirals_pre_cell{kk,1};
    [lia1,locb1] = ismember(spirals_pre_temp(:,1:2),indexSSp_left,'rows');
    spirals_pre_temp2 = spirals_pre_temp(lia1,:);
    ratio_pre2(kk,1) = sum(spirals_pre_temp2(:,4) == -1)./size(spirals_pre_temp2,1);
end
%%
mouseN = 5;
h5f = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
% ratio_pre2 = 1-ratio_pre2;
% ratio_post2 = 1-ratio_post2;
subplot(1,2,1);
scatter(ones(1,mouseN),ratio_pre2,4,'k');
hold on;
scatter(ones(1,mouseN)*2,ratio_post2,4,'k');
hold on;
plot([ones(1,mouseN)*1;ones(1,mouseN)*2],[ratio_pre2';ratio_post2'],'k');
hold on;
yline(0.5,'k--');
ylim([0.4,0.8]);
xticks([1,2]);
xticklabels({'Pre','Post'});
ylabel('CCW Ratio')
[h3,p3] = ttest(ratio_pre2',ratio_post2');
text(1.2,.7,['p = ' num2str(p3)]);

subplot(1,2,2);
scatter(ones(1,mouseN),ratio_pre,4,'k');
hold on;
scatter(ones(1,mouseN)*2,ratio_post,4,'k');
hold on;
plot([ones(1,mouseN)*1;ones(1,mouseN)*2],[ratio_pre';ratio_post'],'k');
yline(0.5,'k--');
ylim([0.4,0.8]);
xticks([1,2]);
xticklabels({'Pre','Post'});
ylabel('CW Ratio')
[h3,p3] = ttest(ratio_pre',ratio_post');
text(1.2,.7,['p = ' num2str(p3)]);

contra_pre_mean = mean(ratio_pre);
contra_post_mean = mean(ratio_post);
contra_pre_sem = std(ratio_pre)./sqrt(5);
contra_post_sem = std(ratio_post)./sqrt(5);
ipsi_pre_mean = mean(ratio_pre2);
ipsi_post_mean = mean(ratio_post2);
ipsi_pre_sem = std(ratio_pre2)./sqrt(5);
ipsi_post_sem = std(ratio_post2)./sqrt(5);
print(h5f, fullfile(save_folder,'Fig5f_SpiralsRatio.pdf'), '-dpdf', '-bestfit', '-painters');
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
hist_bin = 40;
h5d = figure('Renderer', 'painters', 'Position', [100 100 800 600]);
lineColor = 'k'; lineColor1 = 'w';
scale3 = 5/1;
trialN = 4000;
hemi = [];
ax1 = subplot(1,2,1);
[unique_spirals1] = density_color_plot2(spirals_pre_all,hist_bin);        
unique_spirals_unit1 = unique_spirals1(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit1 = unique_spirals_unit1./trialN/5*35;                   % spirals/(mm^2*s)
scatter(unique_spirals1(:,1),unique_spirals1(:,2),...
    3,unique_spirals_unit1,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse');
colormap(ax1,hot);
axis off; axis image;
xlim(ax1,[0,1140]);
ylim(ax1,[0,1320]);
caxis([0,3]);
colorbar;

ax2 = subplot(1,2,2);
[unique_spirals2] = density_color_plot2(spirals_post_all,hist_bin);        
unique_spirals_unit2 = unique_spirals2(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit2 = unique_spirals_unit2./trialN/5*35;                  % spirals/(mm^2*s)
scatter(unique_spirals2(:,1),unique_spirals2(:,2),...
    3,unique_spirals_unit2,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse');
colormap(ax2,hot);
axis off; 
axis image;
xlim(ax2,[0,1140]);
ylim(ax2,[0,1320]);
caxis([0,3]);
colorbar;
print(h5d, fullfile(save_folder,'Fig5d_SpiralsPrePost.pdf'), '-dpdf', '-bestfit', '-painters');