function [hr4e,hr4f] = plotPowerSpectrumParams(T,data_folder,save_folder)
%% plot all power spectrum
data_folder1 =  fullfile(data_folder, 'spirals','spirals_power_spectrum2',...
    'example_traces_005_8Hz');
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
nameList2 = {'VISp','RSP','SSp','MOs'};
freq_value = [0.05,0.1,0.2,0.5,1,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 350;
freq2 = 201; % index of 10Hz
% FR 8Hz at sample 161;
%%
% color1 = cbrewer2('seq','BuPu',15);
color1 = flipud(gray(15));
hr4e = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
subplot(1,1,1);
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder1,[fname '_fft.mat']));
    %%
    psdx_SSp = mean(psdx_mean(:,3:7),2);
    psdx_mean2 = cat(2,psdx_mean(:,1:2),psdx_SSp,psdx_mean(:,8));
    % freq = 0.5 is at index 11
    x = log10(freq1(11:freq2));
    y = log10(psdx_mean2(11:freq2,3));
    y1 = interp1([x(1),x(end)],[y(1),y(end)],x);
    y_ratio = y'-y1;
    plot(x,y_ratio,'color',color1(kk,:));
    hold on;
    y_ratio_all(kk,:) = y_ratio;
    %%
    slope(kk) = (y(end)-y(1))./(x(end)-x(1));
    offset(kk) = y(1);
end
y_ratio_mean = mean(y_ratio_all,1);
y_ratio_sem = std(y_ratio_all,[],1)/sqrt(15);
%
shadedErrorBar(x,y_ratio_mean,y_ratio_sem, 'lineprops', '-r');
freq_value = [0.5,2,4,6,8,10];
log_freq_value = log10(freq_value);
xticks(log_freq_value);      
xticklabels({'0.5','2','4','6','8','10'});
xline(log_freq_value(3),'k--');
xlim([log10(0.5),log10(10)]);
%%
print(hr4e, fullfile(save_folder,'FigR4e_power_diff'),'-dpdf', '-bestfit', '-painters');
%%
freqb = freq1(11:freq2);
peak_bump = y_ratio_all(:,71); % 71 is index for fre1 = 4
% load spiral wave peak density
load(fullfile(data_folder,'spirals\spirals_density',...
    'spiralDensityLinePerSession.mat'))
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
%%
mdl_all = fitlm([peak_bump, slope', offset'],max_density');
%%
mdla = fitlm(max_density,peak_bump);
r2a = mdla.Rsquared.Adjusted;
mdlb = fitlm(max_density,slope);
r2b = mdlb.Rsquared.Adjusted;
mdlc = fitlm(max_density,offset);
r2c = mdlc.Rsquared.Adjusted;
%%
hr4f = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
subplot(1,1,1);
% scatter(max_density,peak_bump);
plot(mdla);
ylim([-0.2,1.2]);
% xlabel('Spiral density (Spirals/mm2*s)');
% ylabel('Spectrum bump amplitude');
xlabel('');
xlabel('');
xticks([0:5]);      
xticklabels(num2cell([0:5]));
print(hr4f, fullfile(save_folder,'FigR4f_spiral_density_vs_bump_amp'),'-dpdf', '-bestfit', '-painters');
%%

