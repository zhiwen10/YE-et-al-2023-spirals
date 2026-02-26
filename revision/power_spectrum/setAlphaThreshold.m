%% manually set and check threshold for each session
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
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    dV = dV(1:50,:);
    % registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
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
    traceAmp3 = traceAmp2(3,:);
    traceFilt3 = traceFilt(3,:);
    figure;
    plot(t,rawTraceV(3,:),'b');
    hold on;
    plot(t,traceFilt3,'k');
    hold on;
    plot(t,traceAmp3,'r');
    %%    
    threshold = 0.006;
    [alpha_epoch2,alpha_binary] = getAlphaEpoch(traceAmp3,threshold);
    traceFilt3b = traceFilt3;
    traceFilt3b(not(alpha_binary)) = nan;
    alpha_ratio = sum(alpha_binary)./numel(alpha_binary);
    %
    figure;
    plot(t,traceFilt3,'k');
    hold on;
    plot(t,traceAmp3,'r');
    hold on;
    yline(threshold,'g--','lineWidth',2);
    hold on;
    plot(t,traceFilt3b,'b');
    %%
    save([fname '_alpha_threshold.mat'],'threshold',...
        'alpha_epoch2','alpha_binary','alpha_ratio');
end