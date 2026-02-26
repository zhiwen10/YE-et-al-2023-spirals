function [f_pxy,Cxy_all,Pxy_all] = coherence2traces(T,data_folder,pixel,behavior)
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder
    fname = [mn '_' tdb '_' num2str(en)];
    subfolder = [mn '_' tdb '_' num2str(en)];
    serverRoot = expPath(mn, td, en);
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    %% registration
    load(fullfile(data_folder,'spirals','rf_tform',...
        [fname '_tform.mat']));  
    Utransformed = imwarp(U,tform,'OutputView',imref2d([1320,1140]));
    trace = squeeze(Utransformed(pixel(1,1),pixel(1,2),1:50))'*V(1:50,:);
    %% load  behavior
    if strcmp(behavior,'pupil')
        tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
        pupil_mean = readNPY(tFile1); 
        % pupilsize = [tUp,pupil_mean];
        pupil_mean2 = pupil_mean(1:2:end);
        if numel(pupil_mean2)<size(t,1)
            pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
        elseif numel(pupil_mean2)>size(t,1)
            pupil_mean2 = pupil_mean2(1:numel(t));
        end 
        signal = pupil_mean2;
    elseif strcmp(behavior,'face')
        load(fullfile(data_folder,'spirals','spirals_index',...
            [fname '_motion_energy.mat']));
        image_energy2(isnan(image_energy2)) = 0;
        signal = image_energy2;
    end
    %% Parameters
    fs = 35;                    % Sampling frequency (Hz)      
    x = signal;  
    y = double(trace)';  
    sizeN = min(numel(x),numel(y));
    x = x(1:sizeN);
    y = y(1:sizeN);
    %% Method 1: Standard mscohere with optimized parameters
    % Window parameters for good 0.1-8 Hz resolution
    window_length = 2^12;       % 4096 samples = ~117 seconds
                               % Gives freq resolution of 35/4096 = 0.0085 Hz
    overlap = window_length/2;  % 50% overlap
    nfft = window_length;      % Same as window length
    % Hamming window for reduced spectral leakage
    window = hamming(window_length);
    % Compute coherence
    [Cxy, f] = mscohere(x, y, window, overlap, nfft, fs);
    %%
    freq_range = (f >= 0.05) & (f <= 8);
    f_roi = f(freq_range);
    Cxy_roi = Cxy(freq_range);
    % Compute cross-spectral density and phase
    [Pxy, f_pxy] = cpsd(x, y, window, overlap, nfft, fs);
    phase_xy = angle(Pxy);
    %%
    Cxy_all(:,kk) = Cxy;
    Pxy_all(:,kk) = Pxy;
end
