function plotSpiralsSymmetryRatioLag(T, data_folder, save_folder)
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
hs = figure('Renderer', 'painters', 'Position', [100 100 900 300]);
subplot(1,3,1);
ms_index = 1;
ratio_all1 = squeeze(ratio_all(ms_index,:,:));
ratio_all1_mean = mean(ratio_all1,1);
ratio_all1_sem = std(ratio_all1,[],1)./sqrt(15);
for i = 1:numel(lags)
    sc = scatter(ones(15,1)*lags(i)/35,ratio_all1(:,i),6,'k','filled');
    hold on;
    er = errorbar(lags/35,ratio_all1_mean,ratio_all1_sem,ratio_all1_sem,'lineWidth',1,'CapSize',8);   
    er.Color = 'g';                            
    hold on;
end 
yline(0.5,'k--');
[h1,p1] = ttest(ratio_all1-0.5);
hold on;
indx1 = find(h1==1,1,'first');
indx2 = find(h1==1,1,'last');
plot([lags(indx1)/35,lags(indx2)/35],[1,1]);
text(-0.2,0.97,[num2str(round(lags(indx1)/35*100)/100), ' to' ,num2str(round(lags(indx2)/35*100)/100)  's']);
title('Left MOp-SSp');
xlabel('Time lag (s)');
ylabel('Proportion of mirrored waves');
xticks([-0.5:0.25:0.5]);
xticklabels([-0.5:0.25:0.5]);
ylim([0.3,1]);

subplot(1,3,2);
ms_index = 3;
ratio_all1 = squeeze(ratio_all(ms_index,:,:));
ratio_all1_mean = mean(ratio_all1,1);
ratio_all1_sem = std(ratio_all1,[],1)./sqrt(15);
for i = 1:numel(lags)
    sc = scatter(ones(15,1)*lags(i)/35,ratio_all1(:,i),6,'k','filled');
    hold on;
    er = errorbar(lags/35,ratio_all1_mean,ratio_all1_sem,ratio_all1_sem,'lineWidth',1,'CapSize',8);   
    er.Color = 'r';                            
    hold on;
end  
yline(0.5,'k--');
[h2,p2] = ttest(ratio_all1-0.5);
hold on;
indx1 = find(h2==1,1,'first');
indx2 = find(h2==1,1,'last');
plot([lags(indx1)/35,lags(indx2)/35],[1,1]);
text(-0.2,0.97,[num2str(round(lags(indx1)/35*100)/100), ' to' ,num2str(round(lags(indx2)/35*100)/100) 's']);
title('SSp left-right');
xlabel('Time lag (s)');
ylabel('Proportion of mirrored waves');
xticks([-0.5:0.25:0.5]);
xticklabels([-0.5:0.25:0.5]);
ylim([0.3,1]);

subplot(1,3,3);
ms_index = 8;
ratio_all1 = squeeze(ratio_all(ms_index,:,:));
ratio_all1_mean = mean(ratio_all1,1);
ratio_all1_sem = std(ratio_all1,[],1)./sqrt(15);
for i = 1:numel(lags)
    sc = scatter(ones(15,1)*lags(i)/35,ratio_all1(:,i),6,'k','filled');
    hold on;
    er = errorbar(lags/35,ratio_all1_mean,ratio_all1_sem,ratio_all1_sem,'lineWidth',1,'CapSize',8);   
    er.Color = 'b';                            
    hold on;
end  
yline(0.5,'k--');
[h3,p3] = ttest(ratio_all1-0.5);
hold on;
indx1 = find(h3==1,1,'first');
indx2 = find(h3==1,1,'last');
plot([lags(indx1)/35,lags(indx2)/35],[1,1]);
text(-0.2,0.97,[num2str(round(lags(indx1)/35*100)/100), ' to' ,num2str(round(lags(indx2)/35*100)/100) 's']);
title('SSp left-left');
xlabel('Time lag (s)');
ylabel('Proportion of mirrored waves');
xticks([-0.5:0.25:0.5]);
xticklabels([-0.5:0.25:0.5]);
ylim([0.3,1]);
%%
print(hs, 'sprial_direction_correlation.pdf','-dpdf', '-bestfit', '-painters');