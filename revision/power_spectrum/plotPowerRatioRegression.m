function hr5beh = plotPowerRatioRegression(T, data_folder,save_folder)
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
hr5beh = figure('Renderer', 'painters', 'Position', [100 100 1200 300]);
ax1 = subplot(1,5,1);
mdla = fitlm(max_density,power_ratio(3,:));
r2a = mdla.Rsquared.Adjusted;
plot(mdla);
ylim([0,50])
xlabel('Spiral density (spirals/mm2*s)');
ylabel('2-8 Hz power ratio');
xlabel('');
xticks([0:5]);      
xticklabels(num2cell([0:5]));
legend('off');

ax2 = subplot(1,5,2);
mdlb = fitlm(max_density,alpha_ratio_all*100');
plot(mdlb);
ylim([0,50]);
ylabel('2-8 Hz epoch ratio');
xlabel('Spiral density (spirals/mm2*s)');
xticks([0:5]);      
xticklabels(num2cell([0:5]));
legend('off');

ax3 = subplot(1,5,3);
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

ax4 = subplot(1,5,4);
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

ax5 = subplot(1,5,5);
peak_bump = y_ratio_all(:,71); % 71 is index for fre1 = 4
mdlc = fitlm(max_density,peak_bump);
plot(mdlc);
ylim([-0.2,1.2]);
xlabel('Spiral density (spirals/mm2*s)');
ylabel('log10(power) (dF/F2)');
xticks([0:5]);      
xticklabels(num2cell([0:5]));
legend('off');
%%
print(hr5beh, fullfile(save_folder,'FigR5beh_spiral_density_vs_power'),'-dpdf', '-bestfit', '-painters');
