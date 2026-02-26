githubdir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals'))); 
addpath(genpath(fullfile(githubdir, 'Pipelines')));     
%% load widefield session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
save_folder = 'E:\spiral_data_share\data\revision2\prediction_motion_energy';
%%
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new2.csv'));
T = T(logical(T.facevideo),:);
T([14,19],:) = [];
session_total = size(T,1);
%%
for kk = 1:size(T,1)
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    serverRoot = expPath(mn, td, en);
    %%
    ops = get_session_info2(T,kk,data_folder);
    expRoot = ops.session_root;
    t = readNPY(fullfile(expRoot,'svdTemporalComponents_corr.timestamps.npy'));
    %%
    image_energy = [];
    image_energy2 = [];    
    %%
    v = VideoReader(fullfile(serverRoot,'face.mp4'));
    numFrames = v.NumFrames;
    batchN = floor(numFrames/1000)+1;
    for i = 1:batchN
        %%
        tic
        frameStart = 1+1000*(i-1);
        frameEnd = 1000+1000*(i-1);
        if i == batchN
            frameEnd = numFrames;
        end
        frames = read(v,[frameStart, frameEnd]);
        images = zeros(size(frames,[1,2,4]));
        image_dff = zeros(size(frames,[1,2,4]));
        for iframe = 1:size(frames,4)
            images(:,:,iframe) = rgb2gray(squeeze(frames(:,:,:,iframe)));
        end
        images = double(images);
        for ii =1:size(images,3)-1
            image_dff(:,:,ii) = images(:,:,ii+1)-images(:,:,ii);
        end
        image_dff(:,:,size(images,3)) = image_dff(:,:,size(images,3)-1);
        image_energy(frameStart:frameEnd) = squeeze(sum(abs(image_dff),[1,2]));
        fprintf('Frame %g/%g; time elapsed %g seconds \n', [frameStart,numFrames, toc])
    end   
    image_energy2 = image_energy(1:2:end);
    if size(image_energy2,2)<size(t,1)
        image_energy2(1,size(image_energy2,2)+1:size(t,1)) = nan;
    elseif size(image_energy2,2)>size(t,1)
        image_energy2 = image_energy2(1:numel(t));
    end       
    image_energy2 = image_energy2';
    %%
    save(fullfile(save_folder,[fname '_motion_energy.mat']),'image_energy2');
end
