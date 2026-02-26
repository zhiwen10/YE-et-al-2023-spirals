function hr5a = plotPowerRatio2(data_folder,save_folder)
%%
load(fullfile(data_folder,'spirals','spirals_power_spectrum2','power_ratio_all_sessions.mat'));
%%
ratio_mean = mean(power_ratio,2);
ratio_sem = std(power_ratio,[],2)./sqrt(15);
%%
color1 = cbrewer2('qual','Set1',9);
%% swap color for VISp and RSP, to match wth example trace
color2 = color1;
color2(2,:) = color1(3,:);
color2(3,:) = color1(2,:);
%%
hr5a = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
index = ones(1,15);
for i = 1:4
    scatter(i*index,power_ratio(i,:),'k');
    hold on;
    errorbar(i,ratio_mean(i),ratio_sem(i),'color',color2(i+1,:),...
        'linewidth',2,'LineStyle','None','Marker','_');
end
hold on;
plot([1*index;2*index],power_ratio(1:2,:),'color',[0.8,0.8,0.8]);
hold on;
plot([2*index;3*index],power_ratio(2:3,:),'color',[0.8,0.8,0.8]);
hold on;
plot([3*index;4*index],power_ratio(3:4,:),'color',[0.8,0.8,0.8]);

xticks([1:4]);
xticklabels(areaNames);
yticks([0:10:60]);
yticklabels(num2cell([0:10:60]));
ylabel('2-8 Hz Power Ratio (%)');
%%
print(hr5a, fullfile(save_folder,'FigR4c_2-8hz_power_ratio'),...
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
load(fullfile(data_folder,'spirals','spirals_density',...
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