function getTaperPowerMap3(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
%% params
psdxAllNorm = zeros(165,143,29,15);
params.tapers=[1,1];
params.Fs=35;
params.fpass=[0.1 8];
params.pad=0;
%%
for kk = 1:size(T,1)
    clearvars -except T formatOut tdAll enAll enexpAll kk ...
        pixSize coords st projectedTemplate2 session_all ...
        sesion_total folder params projectedTemplate1 ...
        data_folder save_folder psdxAllNorm
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
    %% registration
    load(fullfile(data_folder,'spirals','rf_tform',...
        [fname '_tform.mat']));                                            % load atlas transformation matrix tform;
    %%
    clear trialV trialVel1 trialT trialVel
    epoch = 35*2;                                                          % number of timestamps for 2 seconds
    tTotal = numel(t);
    epochN = floor(tTotal/epoch);
    trialV = zeros(200,epochN,epoch);
    for k = 1:epochN
        trialV(:,k,:) = V(:,1+(k-1)*epoch:k*epoch);
        trialT(k,:) = t(1+(k-1)*epoch);
    end
    %% example trace SVD
    nV = size(trialV,1);
    trialN = size(trialV,2);
    trialLength = size(trialV,3);
    thisV = reshape(trialV,size(trialV,1)*size(trialV,2),size(trialV,3));
    thisV = bsxfun(@minus, thisV, mean(thisV, 2)); % mean-subtract each SVD
    [J1,f1]=computePowerMap2(thisV',params);
    nTaper = params.tapers(2);
    J1 = reshape(J1,size(J1,1),nTaper,nV,trialN);
    J1 = permute(J1,[3,1,2,4]);
    J1 = reshape(J1,nV,numel(f1)*nTaper*trialN);
    U1 = U(1:8:end,1:8:end,:)./mimg(1:8:end,1:8:end);
    Ur = reshape(U1,size(U1,1)*size(U1,2),size(U1,3));
    Jpixel = Ur*J1;
    Spixel = conj(Jpixel).*Jpixel;
    Spixel = reshape(Spixel,size(Spixel,1),numel(f1),nTaper*trialN);
    SpixelMean = squeeze(mean(Spixel,3));
    SpixelMean = reshape(SpixelMean,size(U1,1),size(U1,2),[]);
    SpixelMean1 = SpixelMean(:,:,1:end);
    SpixelMean2 = squeeze(mean(mean(SpixelMean,1),2));
    SpixelMean1 = imresize(SpixelMean1,8);
    %%
    psdxMeanTransformed = imwarp(SpixelMean1,tform,...
        'OutputView',imref2d(size(projectedTemplate1))); 
    %%
    psdxMeanTransformed2 = psdxMeanTransformed(1:8:end,1:8:end,:);
    psdxAllNorm(:,:,:,kk) = psdxMeanTransformed2;
%     save(fullfile(save_folder,[fname '_fftSpectrum2']),...
%         'psdxMeanTransformed','freq');
end
%%
freq = f1;
save(fullfile(save_folder,'fftSpectrumAllNew2.mat'),'psdxAllNorm','freq');