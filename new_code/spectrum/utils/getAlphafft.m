function getAlphafft(T,freq,data_folder,save_folder)
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
    clear rawTrace traceFilf traceHilbert tracePhase traceAmp firstFrame ...
        lastFrame traceEpoch traceEpoch1 psdx...
        traceAlpha1 traceNonalpha1 psdx_alpha psdx_nonalpha
    %% session info
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
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));
    %%
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,1:50);
    mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
    %%
    pixel_copy = pixel;
    if kk == 3
        pixel_copy(5,:) = [560,920];
        pixel_copy(6,:) = [655,890];
    end
    pixel_copy = round(pixel_copy/params.downscale);
    %%        
    for i = 1:8
        rawTrace(i,:) = squeeze(Utransformed(pixel_copy(i,1),pixel_copy(i,2),:))'*V(1:50,:);
        rawTrace(i,:) = rawTrace(i,:)./mimgtransformed(pixel_copy(i,1),pixel_copy(i,2));
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
    %% mean psdx acoss epochs and SSp areas
    epochLength = 70;
    [Tall,Ta,Tna,thred1,thred2,traceAlpha,traceNonalpha] = getAmpEpochs3(traceAmp,epochLength,rawTrace);
    traceAlpha = traceAlpha(3:7,:,:);
    traceNonalpha = traceNonalpha(3:7,:,:);
    traceAlpha1 = reshape(traceAlpha,size(traceAlpha,1)*size(traceAlpha,2),[]);
    traceNonalpha1 = reshape(traceNonalpha,size(traceNonalpha,1)*size(traceNonalpha,2),[]);
    [freq1,psdx_alpha,psdx_mean] = fft_spectrum(traceAlpha1);
    [freq1,psdx_nonalpha,non_psdx_mean] = fft_spectrum(traceNonalpha1);
    psdx_alpha_mean = mean(psdx_alpha,2);
    psdx_nonalpha_mean = mean(psdx_nonalpha,2);
    %%
%     h3 = figure('Renderer', 'painters', 'Position', [100 100 250 400]);
%     freq_value = [0.5,2,4,6,8,10];
%     log_freq_value = log10(freq_value);
%     freqN = 21;
%     plot(log10(freq1(2:freqN)),log10(psdx_alpha_mean(2:freqN)),'r');
%     hold on;
%     plot(log10(freq1(2:freqN)),log10(psdx_nonalpha_mean(2:freqN)),'k');
%     xlim([log_freq_value(1),log_freq_value(end)]);
%     xticks(log_freq_value);
%     xticklabels({'0.5','2','4','6','8','10'});
%     xlabel('log10(Frequency)');
%     ylabel('log10(Power) (df/f^2)');
%     print(h3, fullfile(save_folder,[fname '_fft_alpha']),...
%     '-dpdf', '-bestfit', '-painters');
%%
    psdx_alpha_mean_all = cat(2,psdx_alpha_mean_all,psdx_alpha_mean);
    psdx_nonalpha_mean_all = cat(2,psdx_nonalpha_mean_all,psdx_nonalpha_mean);
end
%%
save(fullfile(save_folder,'SSp_fft_all.mat'),...
    'psdx_alpha_mean_all','psdx_nonalpha_mean_all','freq1')
%%
psdx_alpha = mean(psdx_alpha_mean_all,2);
psdx_nonalpha = mean(psdx_nonalpha_mean_all,2);
h3 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
subplot(1,4,1);
freq_value = [0.5,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 21;
plot(log10(freq1(2:freqN)),log10(psdx_alpha(2:freqN)),'r');
hold on;
plot(log10(freq1(2:freqN)),log10(psdx_nonalpha(2:freqN)),'k');
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value);      
xticklabels({'0.5','2','4','6','8','10'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (df/f^2)');
print(h3, fullfile(save_folder,['SSp_fft_mean']),...
    '-dpdf', '-bestfit', '-painters');
