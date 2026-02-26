function hr4c = plotPowerRatio(data_folder,save_folder)
%%
load(fullfile(data_folder,'spirals','spectrum','power_ratio_all_sessions.mat'));
ratio_mean = mean(power_ratio,2);
ratio_sem = std(power_ratio,[],2)./sqrt(15);
%%
hr4c = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
index = ones(1,15);
scatter(1*index,power_ratio(1,:),'k');
hold on;
scatter(2*index,power_ratio(2,:),'k');
hold on;
scatter(3*index,power_ratio(3,:),'k');
hold on;
scatter(4*index,power_ratio(4,:),'k');
hold on;
plot([1*index;2*index],power_ratio(1:2,:),'k');
hold on;
plot([2*index;3*index],power_ratio(2:3,:),'k');
hold on;
plot([3*index;4*index],power_ratio(3:4,:),'k');
hold on;
bar([1,2,3,4],ratio_mean,'edgeColor','k','faceColor','None','linewidth',2);
hold on;
errorbar([1,2,3,4],ratio_mean,ratio_sem,'k','linewidth',2,'LineStyle','none');
xticks([1:4]);
xticklabels(areaNames);
yticks([0:10:60]);
yticklabels(num2cell([0:10:60]));
ylabel('2-8 Hz Power Ratio (%)');
%%
print(hr4c, fullfile(save_folder,'FigR4c_2-8hz_power_ratio'),...
    '-dpdf', '-bestfit', '-painters');
%% 
% compare RSP to VISp, SSP, MOs
count1  =1;
for i = [1,3,4]
    [h1(count1),p1(count1)] = ttest(power_ratio(2,:),power_ratio(i,:));
    count1 = count1+1;
end
% compare SSp to VISp, RSP, MOs
count1  =1;
for i = [1,2,4]
    [h2(count1),p2(count1)] = ttest(power_ratio(3,:),power_ratio(i,:));
    count1 = count1+1;
end    
%% power ratio vs sprial density regression
load(fullfile(data_folder,'spirals\spirals_density',...
    'spiralDensityLinePerSession.mat'))
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
%%
[r1,p1] = corrcoef(max_density, power_ratio(3,:));
mdla = fitlm(max_density,power_ratio(3,:));
r2a = mdla.Rsquared.Adjusted;
% figure;
% subplot(1,1,1);
% plot(mdla);