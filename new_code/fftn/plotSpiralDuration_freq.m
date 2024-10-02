function h1d = plotSpiralDuration_freq(ax,T,freq,data_folder)
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
%% compile duration ratio
for kk = 1:size(T,1)
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en) '.mat'];
    filename = fullfile(data_folder,'spirals\spirals_freq',...
        'spirals_duration',freq_folder,fname);
    load(filename);
    %%
    N_ratio_scramble_all(:,kk) = mean(N_ratio_scramble,2);
    N_ratio_all(:,kk) = N_ratio;
end
%% plot mean and STD
N_ratio_all1 = (N_ratio_all);
N_ratio_scramble_all1 = (N_ratio_scramble_all);
mean_N = mean(N_ratio_all1,2);
std_N = std(N_ratio_all1,[],2);
mean_N_scramble = mean(N_ratio_scramble_all1,2);
std_N_scramble = std(N_ratio_scramble_all1,[],2);

x = [1:50]/35;
errorbar(ax,x,mean_N, std_N,'r',"CapSize",15);
hold on;
errorbar(ax,x,mean_N_scramble,std_N_scramble,'k',"CapSize",15);
xlim([0,15/35]);
ylim([0,1]);
xticks([0, 0.2, 0.4, 0.6,0.8 1.0])
t1 = [0, 0.2, 0.4, 0.6,0.8 1.0];
t2 = cellstr(string(t1));
xticklabels(t2);
xlabel('Spiral duration (ms)');
ylabel('Spiral ratio');
title(freq_folder,'interpreter','none');