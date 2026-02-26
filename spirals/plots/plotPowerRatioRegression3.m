function hs1f = plotPowerRatioRegression2(T, data_folder,save_folder)
%%
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
nameList2 = {'VISp','RSP','SSp','MOs'};
freq_value = [0.05,0.1,0.2,0.5,1,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 350;
freq2 = 201; % index of 10Hz
%%
% load power ratio
load(fullfile(data_folder,'spirals','spirals_power_spectrum2','power_ratio_all_sessions.mat'));
% power ratio vs sprial density regression
load(fullfile(data_folder,'spirals','spirals_density',...
    'spiralDensityLinePerSession.mat'))
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
% [r1,p1] = corrcoef(max_density, power_ratio(3,:));
%%
data_folder1 =  fullfile(data_folder, 'spirals','spirals_power_spectrum2',...
    'example_traces_005_8Hz');
alpha_ratio_all = [];
for kk = 1:size(T,1)
    % session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_power_spectrum2','alpha_threshold',...
        [fname '_alpha_threshold.mat']));
    alpha_ratio_all(kk,1) = alpha_ratio;
    %%
    load(fullfile(data_folder1,[fname '_fft.mat']));
    psdx_SSp = mean(psdx_mean(:,3:7),2);
    psdx_mean2 = cat(2,psdx_mean(:,1:2),psdx_SSp,psdx_mean(:,8));
    % freq = 0.5 is at index 11
    x = log10(freq1(11:freq2));
    y = log10(psdx_mean2(11:freq2,3));
    y1 = interp1([x(1),x(end)],[y(1),y(end)],x);
    y_ratio = y'-y1;
    y_ratio_all(kk,:) = y_ratio;
    slope(kk) = (y(end)-y(1))./(x(end)-x(1));
    offset(kk) = y(1);
end
%%
max_density1 = max_density(1:11);
max_density2 = max_density(12:15);
mean_density1 = mean(max_density1);
mean_density2 = mean(max_density2);
sem_density1 = std(max_density1)./sqrt(11);
sem_density2 = std(max_density2)./sqrt(4);
%%
power_ratio_ssp = power_ratio(3,:);
mean_ratio1 = mean(power_ratio_ssp(1:11));
mean_ratio2 = mean(power_ratio_ssp(12:15));
sem_ratio1 = std(power_ratio_ssp(1:11))./sqrt(11);
sem_ratio2 = std(power_ratio_ssp(12:15))./sqrt(4);
%%
alpha_ratio_all = alpha_ratio_all*100;
mean_alpha_ratio1 = mean(alpha_ratio_all(1:11));
mean_alpha_ratio2 = mean(alpha_ratio_all(12:15));
sem_alpha_ratio1 = std(alpha_ratio_all(1:11))./sqrt(11);
sem_alpha_ratio2 = std(alpha_ratio_all(12:15))./sqrt(4);
%%
mean_alpha_ratio_all = mean(alpha_ratio_all);
sem_alpha_ratio_all = std(alpha_ratio_all)./sqrt(15);

%%
[h0,p0] = ttest2(max_density(1:11),max_density(12:15));
[h1,p1] = ttest2(power_ratio_ssp(1:11),power_ratio_ssp(12:15));
[h2,p2] = ttest2(alpha_ratio_all(1:11),alpha_ratio_all(12:15));
%%
hs1f = figure('Renderer', 'painters', 'Position', [100 100 1200 300]);
ax1 = subplot(1,4,1);
scatter(max_density(1:11),power_ratio_ssp(1:11),'k','filled');
hold on;
errorbar(mean_density1,mean_ratio1,sem_ratio1,sem_ratio1,sem_density1,sem_density1,...
'marker','o','color','k');
hold on;
scatter(max_density(12:15),power_ratio_ssp(12:15),'r','filled');
hold on;
errorbar(mean_density2,mean_ratio2,sem_ratio2,sem_ratio2,sem_density2,sem_density2,...
'marker','o','color','r');
xlim([0,5]);
ylim([0,50]);
xticks([0:5]);      
xticklabels(num2cell([0:5]));
xlabel('Spiral density (spirals/mm2*s)');
ylabel('2-8 Hz power ratio (%)');

ax2 = subplot(1,4,2);
scatter(max_density(1:11),alpha_ratio_all(1:11),'k','filled');
hold on;
errorbar(mean_density1,mean_alpha_ratio1,sem_alpha_ratio1,sem_alpha_ratio1,sem_density1,sem_density1,...
'marker','o','color','k');
hold on;
scatter(max_density(12:15),alpha_ratio_all(12:15),'r','filled');
hold on;
errorbar(mean_density2,mean_alpha_ratio2,sem_alpha_ratio2,sem_alpha_ratio2,sem_density2,sem_density2,...
'marker','o','color','r');
xlim([0,5]);
ylim([0,50]);
ylabel('2-8 Hz epoch ratio (%)');
xlabel('Spiral density (spirals/mm2*s)');
xticks([0:5]);      
xticklabels(num2cell([0:5]));

ax3 = subplot(1,4,3);
kk = 12;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(td,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','spirals_power_spectrum2','alpha_threshold',...
    [fname '_alpha_threshold.mat']));
alpha_ratio_all(kk,1) = alpha_ratio;
load(fullfile(data_folder1,[fname '_fft.mat']));
psdx_SSp = mean(psdx_mean(:,3:7),2);
psdx_mean2 = cat(2,psdx_mean(:,1:2),psdx_SSp,psdx_mean(:,8));
% freq = 0.5 is at index 11
x = log10(freq1(11:freq2));
y = log10(psdx_mean2(11:freq2,3));
y1 = interp1([x(1),x(end)],[y(1),y(end)],x);

plot(x,y,'color','k');
hold on;
plot(x,y1,'color',[0.2,0.2,0.2]);
freq_value = [0.5,2,4,6,8,10];
log_freq_value = log10(freq_value);
xticks(log_freq_value);      
xticklabels({'0.5','2','4','6','8','10'});
% ylim([-0.2,1.2]);
xline(log_freq_value(3),'k--');
xlim([log10(0.5),log10(10)]);
xlabel('log10(frequency)');
ylabel('log10(power) (dF/F2)');

ax4 = subplot(1,4,4);
color1 = flipud(gray(15));
color2 = cbrewer2('seq','OrRd',6);
color3 = [color1(1:11,:);color2(3:6,:)];
for kk = 1:15
    plot(x,y_ratio_all(kk,:),'color',color3(kk,:));
    hold on;
end
freq_value = [0.5,2,4,6,8,10];
log_freq_value = log10(freq_value);
xticks(log_freq_value);      
xticklabels({'0.5','2','4','6','8','10'});
ylim([-0.2,1.2]);
xline(log_freq_value(3),'k--');
xlim([log10(0.5),log10(10)]);
xlabel('log10(frequency)');
ylabel('log10(power) (dF/F2)');
%%
print(hs1f, fullfile(save_folder,'FigS1f_spiral_density_vs_power2'),'-dpdf', '-bestfit', '-painters');
