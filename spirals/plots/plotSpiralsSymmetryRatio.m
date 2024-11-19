function hs7ef = plotSpiralsSymmetryRatio(T, data_folder, save_folder)
%%
load(fullfile(data_folder,'spirals','spirals_symmetry','roiSelection.mat'));
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
    ratio = getSymmetryRatio(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spirals2);
    ratio_all(:,kk) = ratio;
end
%%
mn1 = T.MouseID(:);
area_name = {'M1-S1-left','M1-S1-right','S1-left-right','M1-left-right','M2-left-right','M1-M2-left','M1-M2-right'};
T1 = array2table(ratio_all,'RowNames',area_name,'VariableNames',mn1);
T1.mean = mean(ratio_all,2);
T1.std = std(ratio_all,[],2);
r = normrnd(0,0.05,[1,size(T,1)]);
%%
hs7ef = figure; 
bar(1:7,T1.mean,'faceColor',[0.5, 0.5, 0.5])
for i = 1:7
    hold on;
    scatter(ones(size(T,1),1)'*i+r,T1{i,1:15},'filled');
end
hold on
er = errorbar(1:7,T1.mean,T1.std,T1.std,'lineWidth',2,'CapSize',18);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

xticks([1 2 3 4 5 6 7]);
xticklabels(area_name);
ylabel('matching ratio');
%%
print(hs7ef,fullfile(save_folder, 'Figs7ef_spirals_symmetry_ratio.pdf'),...
    '-dpdf', '-bestfit', '-painters');