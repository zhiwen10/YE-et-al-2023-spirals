function hr3b = plotExamplePowerSpectrum(T,data_folder,save_folder)
%% plot all power spectrum
data_folder1 =  fullfile(data_folder, 'spirals\spectrum\example_traces_005_8Hz');
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
nameList2 = {'VISp','RSP','SSp','MOs'};
freq_value = [0.05,0.1,0.2,0.5,1,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 350;
% FR 8Hz at sample 161;
%%
psdx_all = [];
% color1 = cbrewer2('seq','BuPu',15);
color1 = flipud(gray(15));
hr3b = figure('Renderer', 'painters', 'Position', [100 100 900 230]);
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
        xticks(log_freq_value(1:end-1));
        xticklabels({'0.05','0.1','0.2','0.5','1','2','4','6','8'});
        xlabel('log10(Frequency)');
        ylabel('log10(Power) (df/f^2)');
        ylim([-9,-2]);   
        title(nameList2(i),'Interpreter','None');
    end
end
print(hr3b, fullfile(save_folder,'FigR3b_example_fft_all'),'-dpdf', '-bestfit', '-painters');
