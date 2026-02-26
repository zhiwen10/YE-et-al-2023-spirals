function hr4d = plotAlphaSpectrum(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
%%
pixel(1,:) = [845,835]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,950]; % SSp-m
pixel(6,:) = [550,950]; % SSp-n
pixel(7,:) = [682,905]; % SSp-bfd
pixel(8,:) = [290,700]; % MOs
color1 = cbrewer2('qual','Set3',8);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
%%
psdx_alpha_mean_all = [];
psdx_nonalpha_mean_all = [];
for kk = 1:size(T,1)
    %%
    clear rawTrace rawTraceV traceFilf traceHilbert tracePhase traceAmp firstFrame ...
        lastFrame traceEpoch traceEpoch1 psdx...
        traceAlpha1 traceNonalpha1 psdx_alpha psdx_nonalpha
    % session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals\svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    dV = dV(1:50,:);
    % registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));
    %
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,1:50);
    mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
    %
    pixel_copy = pixel;
    if kk == 3
        pixel_copy(5,:) = [560,920];
        pixel_copy(6,:) = [655,890];
    end
    pixel_copy = round(pixel_copy/params.downscale);
    %%
    load(fullfile(data_folder,'spirals','spirals_power_spectrum2','alpha_threshold',...
        [fname '_alpha_threshold.mat']));
    %%
    alpha_ratio_all(kk,1) = alpha_ratio;
    %%
    mean_all(kk,1) = mean(alpha_epoch2(:,2)-alpha_epoch2(:,1));
    median_all(kk,1) = median(alpha_epoch2(:,2)-alpha_epoch2(:,1));
    %%    
    for i = 1:8
        rawTrace(i,:) = squeeze(Utransformed(pixel_copy(i,1),pixel_copy(i,2),:))'*dV(1:50,:);
        rawTrace(i,:) = rawTrace(i,:)./mimgtransformed(pixel_copy(i,1),pixel_copy(i,2));
    end
    %
     for i = 1:8
        rawTraceV(i,:) = squeeze(Utransformed(pixel_copy(i,1),pixel_copy(i,2),:))'*V(1:50,:);
        rawTraceV(i,:) = rawTraceV(i,:)./mimgtransformed(pixel_copy(i,1),pixel_copy(i,2));
    end
    %%
    Fs = 35;
    rawTrace = double(rawTrace);
    rawTrace = rawTrace -mean(rawTrace ,2); 
    freq = [2,8];
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    traceFilt = filtfilt(f1,f2,rawTrace');
    traceFilt = traceFilt';
    traceHilbert =hilbert(traceFilt);
    tracePhase = angle(traceHilbert);
    traceAmp = abs(traceHilbert);    
    traceAmp2 = movmean(traceAmp,17,2); % smooth mean of 0.5s
    %%
    clear nonalpha_epoch alpha_epoch3 nonalpha_epoch3 raw_trace1 nonalpha_raw_trace1 
    %%
    alpha_size = alpha_epoch2(:,2)-alpha_epoch2(:,1);
    alpha_epoch3(:,1) = alpha_epoch2(alpha_size>=70,1);
    alpha_epoch3(:,2) = alpha_epoch3(:,1)+70;
    %
    for i = 1:size(alpha_epoch3,1)
        raw_trace1(i,:,:) = rawTraceV(3:7,alpha_epoch3(i,1):alpha_epoch3(i,2));
    end
    raw_trace1 = reshape(raw_trace1,size(raw_trace1,1)*size(raw_trace1,2),size(raw_trace1,3));
    [freq1,psdx_alpha,psdx_alpha_mean] = fft_spectrum(raw_trace1);
    %%
    nonalpha_epoch(:,1) = alpha_epoch2(1:end-1,2);
    nonalpha_epoch(:,2) = alpha_epoch2(2:end,1);
    nonalpha_size = nonalpha_epoch(:,2)-nonalpha_epoch(:,1);
    nonalpha_epoch3(:,1) = nonalpha_epoch(nonalpha_size>=70,1);
    nonalpha_epoch3(:,2) = nonalpha_epoch3(:,1)+70;
    
    for i = 1:size(nonalpha_epoch3,1)
        nonalpha_raw_trace1(i,:,:) = rawTraceV(3:7,nonalpha_epoch3(i,1):nonalpha_epoch3(i,2));
    end
    nonalpha_raw_trace1 = reshape(nonalpha_raw_trace1,...
        size(nonalpha_raw_trace1,1)*size(nonalpha_raw_trace1,2),size(nonalpha_raw_trace1,3));
    [freq1,psdx_nonalpha,psdx_nonalpha_mean] = fft_spectrum(nonalpha_raw_trace1);
    %%
    psdx_alpha_mean_all = cat(2,psdx_alpha_mean_all,psdx_alpha_mean);
    psdx_nonalpha_mean_all = cat(2,psdx_nonalpha_mean_all,psdx_nonalpha_mean);
end
%%
psdx_alpha = mean(psdx_alpha_mean_all,2);
psdx_nonalpha = mean(psdx_nonalpha_mean_all,2);
psdx_alpha_sem = std(log10(psdx_alpha_mean_all),[],2)./sqrt(15);
psdx_nonalpha_sem = std(log10(psdx_nonalpha_mean_all),[],2)./sqrt(15);
%%
mean_alpha_ratio = mean(alpha_ratio_all);
sem_alpha_ratio = std(alpha_ratio_all)./sqrt(15);
mean_median_duration = mean(median_all);
sem_median_duration = std(median_all)./sqrt(15);
mean_mean_duration = mean(mean_all);
sem_mean_duration = std(mean_all)./sqrt(15);
%%
hr4d = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
subplot(1,3,[1,2]);
freq_value = [0.5,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 21;
shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_alpha(2:freqN)),psdx_alpha_sem(2:freqN), 'lineprops', '-r');
% plot(log10(freq1(2:freqN)),log10(psdx_alpha(2:freqN)),'r');
hold on;
shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_nonalpha(2:freqN)),psdx_nonalpha_sem(2:freqN), 'lineprops', '-k');
% plot(log10(freq1(2:freqN)),log10(psdx_nonalpha(2:freqN)),'k');
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value);      
xticklabels({'0.5','2','4','6','8','10'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (df/f^2)');
subplot(1,3,3);
scatter(ones(15,1),alpha_ratio_all,'k');
hold on;
bar(1,mean_alpha_ratio);
hold on;
errorbar(1,mean_alpha_ratio,sem_alpha_ratio,'linewidth',2,'LineStyle','None','Marker','_');
ylim([0,0.5]);
yticks([0:0.1:0.5]);
yticklabels(num2cell([0:0.1:0.5]));
ylabel('2-8Hz epoch ratio');
%%
print(hr4d, fullfile(save_folder,'FigR4d_alpha_bump'),...
    '-dpdf', '-bestfit', '-painters');