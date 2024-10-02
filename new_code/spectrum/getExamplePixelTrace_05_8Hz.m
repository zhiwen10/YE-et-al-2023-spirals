function getExamplePixelTrace_05_8Hz(T,data_folder,save_folder)
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
psdx_all1 = [];
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
    %%
    epochLength = 70;
    epochs = floor(size(rawTrace,2)./epochLength);
    for i = 1:8
        for epoch_i = 1:epochs
            firstFrame(epoch_i,1) = 1+(epoch_i-1)*epochLength;
            lastFrame(epoch_i,1) = epoch_i*epochLength;
            traceEpoch(i,epoch_i,:) = rawTrace(i,firstFrame(epoch_i):lastFrame(epoch_i));
        end
    end
    %%
    traceEpoch1 = reshape(traceEpoch,size(traceEpoch,1)*size(traceEpoch,2),[]);
    [freq1,psdx,~] = fft_spectrum(traceEpoch1);
    psdx = reshape(psdx,36,size(traceEpoch,1),size(traceEpoch,2));
    psdx_mean = mean(psdx,3);
    %%
    freq_value = [0.5,2,4,6,8];
    log_freq_value = log(freq_value);
    freqN = 17;
    h1 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
    for i = 1:8
        subplot(2,4,i);
        plot(log(freq1(2:freqN)),log(psdx_mean(2:freqN,i)),'r');
        xlim([log_freq_value(1),log_freq_value(end)]);
        xticks(log_freq_value);
        xticklabels({'0.5','2','4','6','8'});
        xlabel('log(Frequency)');
        ylabel('log(Power) (df/f^2)');
        ylim([-18,-6]);
        title(nameList(i),'Interpreter','None');
    end
    print(h1, fullfile(save_folder,[fname '_fft_all']),...
    '-dpdf', '-bestfit', '-painters');
    %%
    epochLength = 70;
    [Tall,Ta,Tna,thred1,thred2,traceAlpha,traceNonalpha] = getAmpEpochs3(traceAmp,epochLength,rawTrace);
    %%
    traceAlpha1 = reshape(traceAlpha,size(traceAlpha,1)*size(traceAlpha,2),[]);
    traceNonalpha1 = reshape(traceNonalpha,size(traceNonalpha,1)*size(traceNonalpha,2),[]);
    [freq1,psdx_alpha,~] = fft_spectrum(traceAlpha1);
    psdx_alpha = reshape(psdx_alpha,36,size(traceAlpha,1),size(traceAlpha,2));
    psdx_alpha_mean = mean(psdx_alpha,3);
    %%
    [freq1,psdx_nonalpha,~] = fft_spectrum(traceNonalpha1);
    psdx_nonalpha = reshape(psdx_nonalpha,36,size(traceNonalpha,1),size(traceNonalpha,2));
    psdx_nonalpha_mean = mean(psdx_nonalpha,3);
    %%
    freq_value = [0.5,2,4,6,8];
    log_freq_value = log(freq_value);
    freqN = 17;
    h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
    for i = 1:8
        subplot(2,4,i);
        plot(log(freq1(2:freqN)),log(psdx_alpha_mean(2:freqN,i)),'r');
        hold on;
        plot(log(freq1(2:freqN)),log(psdx_nonalpha_mean(2:freqN,i)),'k');
        xlim([log_freq_value(1),log_freq_value(end)]);
        xticks(log_freq_value);
        xticklabels({'0.5','2','4','6','8'});
        xlabel('log(Frequency)');
        ylabel('log(Power) (df/f^2)');
        ylim([-18,-6]);
        title(nameList(i),'Interpreter','None');
    end
    print(h2, fullfile(save_folder,[fname '_fft_alpha']),...
    '-dpdf', '-bestfit', '-painters');
    %% 
    close all;
    %%
    save(fullfile(save_folder,[fname '_fft.mat']),...
        'freq1', 'psdx_mean', 'psdx_alpha_mean', 'psdx_nonalpha_mean');
    %%
    psdx_all = cat(3,psdx_mean,psdx_alpha_mean,psdx_nonalpha_mean);
    psdx_all1 = cat(4,psdx_all1,psdx_all);
end
writeNPY(psdx_all1,fullfile(save_folder,'fft.npy'));
writeNPY(freq1,fullfile(save_folder,'frequency.npy'));