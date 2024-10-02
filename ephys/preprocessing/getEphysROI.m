function getEphysROI(T,data_folder,save_folder)
%%
params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
params.padding = 2*params.halfpadding;
roi_exist = zeros(size(T,1),1);
for kk = 1:size(T,1)
    clearvars -except T kk params data_folder save_folder roi_exist count1 
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)]; 
    mimg = readNPY(fullfile(ops.session_root, 'meanImage.npy'));
    %% apply mask, this helps speed up spiral detection later
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    %%
    fname1 = [fname '_roi.mat'];
    filename = fullfile(save_folder,fname1);
    if exist(filename, 'file') == 2
        roi_exist(kk) = 1;
    else
        roi_exist(kk) = 0;
        figure; 
        ax1 = imagesc(mimg2);
        roi = drawpolygon;
        save(fullfile(save_folder,fname1),'roi');
        close all;
    end
end

