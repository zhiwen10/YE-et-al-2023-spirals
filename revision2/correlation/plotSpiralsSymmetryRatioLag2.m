function plotSpiralsSymmetryRatioLag2(T, data_folder, save_folder)
%%
load(fullfile(data_folder,'spirals','spirals_symmetry','roiSelection.mat'));
%%
lags = -17:17;
ratio_all = [];
for kk = 1:15
    %%
    clear spiralsT 
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals2 = cell2mat(groupedCells); 
    %%
    filteredSpirals2(filteredSpirals2(:,4) == -1,4) = 0;
    spirals2 = filteredSpirals2;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,spirals2(:,1),spirals2(:,2));
    spirals2(:,1:2) = spiralsT;        
    for i = 1:numel(lags)
        clear ratio
        ratio = getSymmetryRatioLag(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spirals2,lags(i));
        ratio2 = getSymmetryRatioLag2(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spirals2,lags(i));
        ratio = [ratio;ratio2];
        ratio_all(:,kk,i) = ratio;
    end
end
%%
ratio_all(8,:,:) = 1-ratio_all(8,:,:);
%%
pair_names = {'AP','HEMI','AUTO'};
color1 = {'g','r','b'};
hs = figure('Renderer', 'painters', 'Position', [100 100 900 250]);
ms_all = [1,3,8];
for m = 1:3
    subplot(1,3,m);
    ms_index = ms_all(m);
    ratio_all1 = squeeze(ratio_all(ms_index,:,:));
    ratio_all1_mean = mean(ratio_all1,1);
    ratio_all1_sem = std(ratio_all1,[],1)./sqrt(15);
    for i = 1:15
        plot(lags/35,ratio_all1(i,:),'k');
        hold on;
    end
    shadedErrorBar(lags/35, ratio_all1_mean, ratio_all1_sem, 'lineprops', ['-' color1{m}]);
    yline(0.5,'k--');
    hold on;
    xline(0,'k--');
    [h1,p1] = ttest(ratio_all1-0.5);
    hold on;
    indx1 = find(h1==1,1,'first');
    indx2 = find(h1==1,1,'last');
    plot([lags(indx1)/35,lags(indx2)/35],[1,1]);
    text(-0.2,0.97,[num2str(round(lags(indx1)/35*100)/100), ' to' ,num2str(round(lags(indx2)/35*100)/100)  's']);
    title(pair_names{m});
    xlabel('Time lag (s)');
    ylabel('Proportion of mirrored waves');
    xticks([-0.5:0.25:0.5]);
    xticklabels([-0.5:0.25:0.5]);
    ylim([0.3,1]);
end
print(hs, 'sprial_direction_correlation.pdf','-dpdf', '-bestfit', '-painters');