function hs13g = plotSpiralsMatchingRate_Arousal(T,data_folder,save_folder)
%%
load(fullfile(data_folder,'ephys','spirals_compare',...
    'spiral_compare_sessions_arousal.mat'));
%%
T = T(logical(T.facevideo),:);
T([14,19],:) = [];
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
match_all = [];
match_perm = [];
for kk = 1:size(spiral_left_match_all,1)
    clear match1 match2 match_temp
    left_temp = spiral_left_match_all{kk,1};
    match1 = left_temp(:,11);
    right_temp = spiral_right_match_all{kk,1};
    match2 = right_temp(:,11);
    match_temp = [match1;match2];   
    match_all(kk,1) = numel(match_temp);
    match_all(kk,2) = sum(match_temp)/numel(match_temp);
    
    clear match1p match2p match_tempp
    left_temp = spiral_left_match_all{kk,2};
    match1p = left_temp(:,11);
    right_temp = spiral_right_match_all{kk,2};
    match2p = right_temp(:,11);
    match_tempp = [match1p;match2p];   
    match_perm(kk,1) = numel(match_tempp);
    match_perm(kk,2) = sum(match_tempp)/numel(match_tempp);
end
%%
color2 = {'g','r','c','m'};
ratio_all = {};
ratio_perm = {};
hs13g = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
count1 = 1;
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    clear indx T2
    indx = contains(T.Area,area(1,kk));
    T2 = T(indx,:);
    match_all_temp = match_all(indx,:);
    match_perm_temp = match_perm(indx,:);
    ratio_all{count1} = match_all_temp;
    ratio_perm{count1} = match_perm_temp;
    subplot(1,3,count1);
    scatter(ones(size(match_all_temp(:,2))),match_all_temp(:,2),18,color2{kk},'filled');
    hold on;
    scatter(2*ones(size(match_perm_temp(:,2))),match_perm_temp(:,2),18,'k','filled');
    hold on;
    plot([ones(size(match_all_temp(:,2))),2*ones(size(match_perm_temp(:,2)))]',[match_all_temp(:,2),match_perm_temp(:,2)]','k');
    hold on;
    % yline(0.5,'--k');
    plot([1,2],[0.5,0.5],'k--');
    ylim([0.3,0.9]);
    [ha2(count1),pa2(count1)] = ttest(match_all_temp(:,2),match_perm_temp(:,2));
    mean_all(count1,1) = mean(match_all_temp(:,2));
    sem_all(count1,1) = std(match_all_temp(:,2))./sqrt(size(match_perm_temp,1));
    mean_perm(count1,1) = mean(match_perm_temp(:,2));
    sem_perm(count1,1) = std(match_perm_temp(:,2))./sqrt(size(match_perm_temp,1));
    xlabel('Spiral counts');
    ylabel('Spiral matching ratio');
    text(1.5,0.7,['p =' num2str(pa2(count1))]);
    yticks([0.3:0.1:0.9])
    yticklabels({'0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
    xlim([0.8,2.2]);
    
    count1 = count1+1;
end
%%
print(hs13g,fullfile(save_folder, 'FigS13g_spirals_matching_ratio_arousal2.pdf'),...
    '-dpdf', '-bestfit', '-painters');