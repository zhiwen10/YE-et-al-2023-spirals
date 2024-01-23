%% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath(fullfile(githubdir2, 'allenCCF')));
addpath(genpath(fullfile(githubdir2, 'AP_scripts_cortexlab-master')));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection'))
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%% only select cortex in the atlas
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'ara_tools-master')))
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
%%
template1 =zeros(1320,1140);
for k1 = 1:1320
    for k2 = 1:1140
        secondIndex = find(av(k1,:,k2)>1,60,'first');
        if not(isempty(secondIndex))
            template1(k1,k2) = tv(k1,secondIndex(end),k2);
        end
    end
end
%%
template2 = template1(1:8:end,1:8:end);
%%
clear tv av
%%
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
% point(8,:) = [110 104]; %% VISp
point(8,:) = [115 105]; %% VISp
point(7,:) = [97 79]; %% RSP
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% cortex surface outline
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
%%
session_all = find(T.use);
session_total = numel(session_all);
%%
ffolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\ReducedRankRregression\hemi-2';
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
count = 1;
for kk = 12:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    mimg = readNPY(fullfile(serverRoot, 'blue', ['meanImage.npy']));
    fname = [mn '_' tdb '_' num2str(en) '-hemi.mat'];
    load(fullfile(ffolder,fname));
    %%
    folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
    dfolder = fullfile(folder,mn,td,num2str(en));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    %%
    % Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %% color overlay a cycle
    color1 = hsv(8);
    mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
    iplot = kk-11;
    axx4 = subplot(6,2,iplot)
    im2 = imagesc(mimgtransformed2);
    colormap(gray)
    BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
    set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');

    hold on;
    for kkk = 1:8
        TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
        image_all(:,:,kkk,count) = TheColorImage;
        % overlay kernel with colormap    
        maxI = max(TheColorImage(:));
        % TheGrayscaleImage = squeeze(mimg);
        % hb1 = imshow(TheGrayscaleImage,[]);
        % colormap(ax1,gray);
        imageSize = size(TheColorImage);  
        colorIntensity = TheColorImage/maxI;
        colorIntensity = colorIntensity.^2;
    %     maxIntensity = max(colorIntensity(:));
    %     colorIntensity(colorIntensity<0.6*maxIntensity) = 0;
        thisColor = color1(kkk,:);
        thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
            thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
        hold on
        hb2(kkk) = imshow(thisColorImage);
        hold off

        set(hb2(kkk),'AlphaData',colorIntensity);
        axis image; 
        hold on;
        scatter(axx4,point(kkk,2),point(kkk,1),8,...
            'MarkerFaceColor',color1(kkk,:),...
            'MarkerEdgeColor','None');
        hold on;
    end
    hold on;
    % overlayOutlines(coords,scale2,'w');
    scale3 = 5/8;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);

    set(gca,'Ydir','reverse')
    axis on; axis off
    %%
    title(fname(1:end-7), 'fontsize',6,'Interpreter', 'none');
    %%
    count = count+1;
end
%%
ax2 = subplot(6,2,5);
colormap(ax2,color1);
cb = colorbar('YTickLabel', nameList,'YTick',0.06:0.13:1.05);
cb.TickLabelInterpreter = 'none';
axis off;
originalSize1 = get(ax2, 'Position');
% set(ax2, 'Position', [originalSize1(1)-0.2,originalSize1(1)-0.2,originalSize1(3)*0.6, originalSize1(4)*0.3]);
%%
print(h1, 'hemi-map-g8', '-dpdf', '-bestfit', '-painters');
%%
count = 1;
ffolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\ReducedRankRregression\hemi-2';
for kk = 1:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    mimg = readNPY(fullfile(serverRoot, 'blue', ['meanImage.npy']));
    fname = [mn '_' tdb '_' num2str(en) '-hemi.mat'];
    load(fullfile(ffolder,fname));
    for kkk = 1:8
        TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
        image_all(:,:,kkk,count) = TheColorImage;
    end
    %%
    count = count+1;
end
%%
save('hemi_kernel_map.mat','image_all');
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
axx4 = subplot(1,1,1);
im2 = imagesc(template2);
colormap(axx4,gray);
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
image_all_mean = mean(image_all,4);    
for kkk = 1:8
    TheColorImage = image_all_mean(:,:,kkk);
    % overlay kernel with colormap    
    maxI = max(TheColorImage(:));
    % TheGrayscaleImage = squeeze(mimg);
    % hb1 = imshow(TheGrayscaleImage,[]);
    % colormap(ax1,gray);
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage/maxI;
    colorIntensity = colorIntensity.^2;
%     maxIntensity = max(colorIntensity(:));
%     colorIntensity(colorIntensity<0.6*maxIntensity) = 0;
    thisColor = color1(kkk,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on
    hb2(kkk) = imshow(thisColorImage);
    hold off

    set(hb2(kkk),'AlphaData',colorIntensity);
    axis image; 
    hold on;
    scatter(axx4,point(kkk,2),point(kkk,1),8,...
        'MarkerFaceColor',color1(kkk,:),...
        'MarkerEdgeColor','None');
    hold on;
end
hold on;
% overlayOutlines(coords,scale2,'w');
scale3 = 5/8;
plotOutline(maskPath(1:3),st,atlas1,[],scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
set(gca,'Ydir','reverse')
axis on; axis off
%%
print(h1, 'hemi-map-all-mean', '-dpdf', '-bestfit', '-painters');
%% ap_division
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
axx4 = subplot(1,1,1);
im2 = imagesc(template2);
colormap(axx4,gray);
scale = 8;
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled'); 
scale3 = 5/8;
hold on;
lineColor = 'k';
plotOutline(maskPath(1:3),st,atlas1,[],scale3,lineColor);
plotOutline(maskPath(4:11),st,atlas1,[],scale3,lineColor);
% plotOutline(maskPath(5),st,atlas1,[],scale3);
% plotOutline(maskPath(6:11),st,atlas1,[],scale3);
set(gca,'Ydir','reverse')
axis on; axis off
axis image;
print(h1, 'AP_division', '-dpdf', '-bestfit', '-painters');
%% hemi_division
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
axx4 = subplot(1,1,1);
im2 = imagesc(template2);
colormap(axx4,gray);
scale = 8;
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled'); 
scale3 = 5/8;
hold on;
lineColor = 'k';
plotOutline(maskPath(1:11),st,atlas1,-1,scale3,lineColor);
plotOutline(maskPath(1:11),st,atlas1,1,scale3,lineColor);
% plotOutline(maskPath(5),st,atlas1,[],scale3);
% plotOutline(maskPath(6:11),st,atlas1,[],scale3);
set(gca,'Ydir','reverse')
axis on; axis off
axis image;
print(h1, 'hemi_division', '-dpdf', '-bestfit', '-painters');
%%
image_all2 = squeeze(image_all(:,:,8,:));
maxIntensity = max(image_all2(:));
minIntensity = min(image_all2(:));
figure;
for i = 1:15
    subplot(4,4,i)
    imagesc(squeeze(image_all2(:,:,i)))
    colormap(gray)
    % caxis([minIntensity,maxIntensity])
end
%%
ffolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\ReducedRankRregression\AP-3';
for kk = 1:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    mimg = readNPY(fullfile(serverRoot, 'blue', ['meanImage.npy']));
    fname = [mn '_' tdb '_' num2str(en) '-AP.mat'];
    load(fullfile(ffolder,fname),'explained_var5');
    var_all_ap(:,kk) =  explained_var5;
end
mean_var_ap_all = mean(var_all_ap,2);
std_var_ap_all = std(var_all_ap,0,2);
sem_var_ap_all = std_var_ap_all/size(var_all_ap, 2);
%%
ffolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\ReducedRankRregression\hemi-3';
for kk = 1:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    mimg = readNPY(fullfile(serverRoot, 'blue', ['meanImage.npy']));
    fname = [mn '_' tdb '_' num2str(en) '-hemi.mat'];
    load(fullfile(ffolder,fname),'explained_var5');
    var_all_hemi(:,kk) =  explained_var5;
end
mean_var_hemi_all = mean(var_all_hemi,2);
std_var_hemi_all = std(var_all_hemi,0,2);
sem_var_hemi_all = std_var_hemi_all/size(var_all_hemi, 2);
%%
h1= figure;
subplot(1,2,1);
shadedErrorBar(1:50, mean_var_ap_all, sem_var_ap_all, 'lineprops', '-r')
ylim([0.65,1.02])

subplot(1,2,2);
shadedErrorBar(1:50, mean_var_hemi_all, sem_var_hemi_all, 'lineprops', '-g')
ylim([0.65,1.02])

%%
print(h1, 'prediction_accuracy_rank', '-dpdf', '-bestfit', '-painters');
%%
for kk = 1:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    fname1 = load(fullfile(serverRoot, 'blue', ['dataSummary.mat']));
    cumsumSv = cumsum(Sv);
    Sv_ratio(:,kk) = cumsumSv/totalVar;
end
%%
mean_Sv_ratio = mean(Sv_ratio,2);
std_Sv_ratio = std(Sv_ratio,0,2);
sem_Sv_ratio = std_Sv_ratio/sqrt(size(Sv_ratio, 2));
h1= figure;
shadedErrorBar(1:500, mean_Sv_ratio(1:500), sem_Sv_ratio(1:500), 'lineprops', '-r')
xline(50,'--');
text(100,0.96,num2str(round(mean_Sv_ratio(50)*100)/100));
print(h1, 'data_variance', '-dpdf', '-bestfit', '-painters');