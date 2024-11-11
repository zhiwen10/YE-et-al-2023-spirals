%% plot all power spectrum
data_folder1 =  fullfile(data_folder, 'spirals\spectrum\example_traces_01_8Hz');
%%
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
nameList2 = {'VISp','RSP','SSp','MOs'};
%%
freq_value = [0.05,0.1,0.2,0.5,1,2,4,6,8,16];
log_freq_value = log10(freq_value);
freqN = 350;
% FR 8Hz at sample 161;
%%
Tfooof = readtable(fullfile(data_folder1,'fooof_15_4_7.csv'));
offset = reshape(Tfooof.offset,4,15);
exponent = reshape(Tfooof.exponent,4,15);
cf = reshape(Tfooof.cf_0,4,15);
pw = reshape(Tfooof.pw_0,4,15);
indx = find(isnan(cf));
[row,col] = ind2sub(size(cf),indx);
col1 = unique(col);
offset(:,col1) = [];exponent(:,col1) = [];
cf(:,col1) = []; pw(:,col1) = [];
%%
psdx_all = [];
color1 = cbrewer2('seq','BuPu',15);
h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
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
    %%
    psdx_all = cat(3,psdx_all,psdx_mean2);
    %%
    for i = 1:4
        subplot(1,4,i);
        plot(log10(freq1(2:freqN)),log10(psdx_mean2(2:freqN,i)),'color',color1(kk,:));
        hold on;
        xlim([log_freq_value(1),log_freq_value(end)]);
        xticks(log_freq_value);
        xticklabels({'0.05','0.1','0.2','0.5','1','2','4','6','8','16'});
        xlabel('log10(Frequency)');
        ylabel('log10(Power) (df/f^2)');
        ylim([-9,-2]);   
        title(nameList2(i),'Interpreter','None');
    end
end
%%
print(h2, fullfile(save_folder,'fft_all'),'-dpdf', '-bestfit', '-painters');
close all;
%%
TfooofExample = readtable(fullfile(data_folder1,'fooofExample_4_7.csv'));
h3 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
x = log_freq_value(3:end-1);
% session info
kk = 7;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(td,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder1,[fname '_fft.mat']));
psdx_SSp = mean(psdx_mean(:,3:7),2);
psdx_mean2 = cat(2,psdx_mean(:,1:2),psdx_SSp,psdx_mean(:,8));
for i = 1:4
    subplot(1,4,i);
    plot(log10(freq1(2:freqN)),log10(psdx_mean2(2:freqN,i)),'color',color1(kk,:));
    hold on;
    y = TfooofExample.offset(i)-TfooofExample.exponent(i)*x;
    plot(x,y,'k--');
    xlim([log_freq_value(1),log_freq_value(end)]);
    xticks(log_freq_value);
    xticklabels({'0.05','0.1','0.2','0.5','1','2','4','6','8','16'});
    xlabel('log10(Frequency)');
    ylabel('log10(Power) (df/f^2)');
    ylim([-9,-2]);   
    title(nameList2(i),'Interpreter','None');
end
print(h3, fullfile(save_folder,'ZYE12_example'),'-dpdf', '-bestfit', '-painters');
%%
load(fullfile(data_folder, 'spirals\spectrum','spirals_fft_radius.mat'));
spirals_sum = sum(count_sample_control(:,4:10),2);
spirals_sum(col1) = [];
xticksValue = [0:2:6];
xticksLabel = string(xticksValue);
h4 = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
for areaN = 1:4
    subplot(4,4,areaN);
    scatter(spirals_sum,offset(areaN,:),6,'k','filled');
    xlim([0,6]); ylim([-6,-4]);
    xlabel('Peak spiral density'); ylabel('Offset');
    xticks(xticksValue);
    xticklabels(xticksLabel);
    
    subplot(4,4,areaN+4);
    scatter(spirals_sum,exponent(areaN,:),6,'k','filled');
    xlim([0,6]); ylim([1.0,3]);
    xlabel('Peak spiral density'); ylabel('Exponent');
    xticks(xticksValue);
    xticklabels(xticksLabel);
    
    subplot(4,4,areaN+8);
    scatter(spirals_sum,cf(areaN,:),6,'k','filled');
    xlim([0,6]); ylim([0,6]);
    xlabel('Peak spiral density'); ylabel('Peak freq');
    xticks(xticksValue);
    xticklabels(xticksLabel);
    
    subplot(4,4,areaN+12);
    scatter(spirals_sum,pw(areaN,:),6,'k','filled');
    xlim([0,6]); ylim([0,1]);
    xlabel('Peak spiral density'); ylabel('Peak power');
    xticks(xticksValue);
    xticklabels(xticksLabel);    
end
print(h4, fullfile(save_folder,'fooof_params'),'-dpdf', '-bestfit', '-painters');
% close all;
%%
figure;
scatter(1,cf(1,:),'k');
hold on;
scatter(3,cf(3,:),'k');
hold on;
scatter(4,cf(4,:),'k');
% plot([ones(1,13);ones(1,13)*2],[cf(1,:);cf(2,:)],'k');
% hold on;
plot([ones(1,13)*1;ones(1,13)*3],[cf(1,:);cf(3,:)],'k');
hold on;
plot([ones(1,13)*3;ones(1,13)*4],[cf(3,:);cf(4,:)],'k');
%% %% unique variance explained
% offsetz = zscore(offset(3,:),[],2);
% exponentz = zscore(exponent(3,:),[],2);
% cfz = zscore(cf(3,:),[],2);
% pwz = zscore(pw(3,:),[],2);
% pwz1 = zscore(pw(1,:),[],2);
% regressor = [offsetz;exponentz;cfz;pwz;pwz1]';
regressor = [exponent(3,:);cf(3,:)]';
mdl = fitlm(regressor,spirals_sum)
%%
r2_total = mdl.Rsquared.Ordinary;