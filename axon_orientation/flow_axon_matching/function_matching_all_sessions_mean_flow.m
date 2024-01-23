%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\nrrdread'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%% get atlas mask and outlines
% load coords for atlas
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
% load projectedAtlas and projectedTemplate
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
ssp_bfd_point = [70,95];
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
for kk = 1:15
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %% SVD, plot example trace, overlay with pupil and stim time
    serverRoot = expPath(mn, td, en);
    [U,V,t,mimg] = get_wf_svd1(serverRoot);
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %%
    dfolder = fullfile(folder,mn,td,num2str(en));
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(dfolder,[fname '_spirals_group_fftn.mat']));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    %%
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
    mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
    BW = logical(projectedAtlas1);
    BW1 = BW(1:params.downscale:end,1:params.downscale:end);
    %%
    [row,col] = find(BW);
    brain_index = [col,row];
    %%
    clear spiralsT
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,filteredSpirals(:,1),filteredSpirals(:,2));    
    filteredSpirals(:,1:2) = round(spiralsT); 
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');
    filteredSpirals = filteredSpirals(lia,:);
    %%
    clear filteredSpirals1
    filteredSpirals1 = filteredSpirals;
    % filteredSpirals1 = filteredSpirals(filteredSpirals(:,3)>=70,:);
    filteredSpirals1 = filteredSpirals1(filteredSpirals1(:,4)==1,:);
    location_index = (filteredSpirals1(:,1)>=800&filteredSpirals1(:,1)<=900&filteredSpirals1(:,2)>=500&filteredSpirals1(:,2)<=650);
    filteredSpirals1 = filteredSpirals1(location_index,:);
    %%
    frame_count = max(filteredSpirals(:,5));
    hist_bin = 40;
    pixSize = 0.01; % mm/pix
    pixArea = pixSize^2;
    %%
    useGPU = 0;
    spiral_n = size(filteredSpirals1,1);
    vxRaw_all = []; vyRaw_all =[];
    vxRaw_all = zeros(165,143,spiral_n); vyRaw_all= zeros(165,143,spiral_n);
    spiral_phase_all = zeros(165,143,2,spiral_n);
    for i = 1:spiral_n
        tracePhase_current  = [];
        frameStart = filteredSpirals1(i,5);
        frameEnd = frameStart+100; 
        frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
        dV1 = dV(:,frameTemp);
        [~,~,tracePhase1] = spiralPhaseMap4(Utransformed,dV1,t,params,rate);
        tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
        tracePhase_current = tracePhase1(:,:,1:2);
        spiral_phase_all(:,:,:,i) = tracePhase_current;        
    end
    %%
    spiral_phase_all_norm = zeros(165,143,2,spiral_n);
    for i = 1:spiral_n
        spiral_phase_temp = squeeze(spiral_phase_all(:,:,:,i));
        spiral_phase_temp1 = spiral_phase_temp-spiral_phase_temp(70,95);
        spiral_phase_temp2 = angle(exp(1i*(spiral_phase_temp1))); 
        spiral_phase_all_norm(:,:,:,i) = spiral_phase_temp2;
    end
    %%
    tracePhase_mean = circ_mean(spiral_phase_all_norm,[],4);
    tracePhase_mean = permute(tracePhase_mean,[3 1 2]);    
    [vxRaw,vyRaw] = HS_flowfield(tracePhase_mean,useGPU);
    %%
    vxRaw_mean = squeeze(vxRaw); vyRaw_mean = squeeze(vyRaw);
    vxRaw2 = nan(size(vxRaw_mean)); vyRaw2 = nan(size(vyRaw_mean));
    skip = 3; zoom_scale = 3;
    vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
    scale2 = 8;
    %%
    h = figure;
    subplot(1,1,1);
    im = imagesc(squeeze(tracePhase_mean(1,:,:)));
    colormap(colorcet('C06'));
    axis image; axis off;
    set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
    hold on;
    imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
    print(h, [fname '_mean_flow2'], '-dpdf', '-bestfit', '-painters');
    save([fname '_mean_flow2'],'spiral_phase_all_norm','spiral_phase_mean','vxRaw_all','vyRaw_all','vxRaw_mean','vyRaw_mean');
    close all;
end