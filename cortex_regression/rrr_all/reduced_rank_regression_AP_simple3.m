%%
% add path for dependencies 
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
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
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
%% right SSp index
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
%% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
%%
session_all = find(T.use);
session_total = numel(session_all);
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
    [U,V,t,mimg] = get_wf_svd1(serverRoot);
    U = U(:,:,1:50);
    V = V(1:50,:);
    % dV = [zeros(size(V,1),1) diff(V,[],2)];
    %%
    folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
    fname = [mn '_' tdb '_' num2str(en)];
    dfolder = fullfile(folder,mn,td,num2str(en));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    %%
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %% mask and Kernel regression map for SSp and MO
    hemi = 'right';
    [indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_SSp,Vnew_SSp] = redoSVD(UselectedSSp,V(:,1:58000));
    USSp = Unew_SSp(:,1:50);
    VSSp = Vnew_SSp(1:50,:);

    [indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_MO,Vnew_MO] = redoSVD(UselectedMO,V(:,1:58000));
    UMO = Unew_MO(:,1:50);
    VMO = Vnew_MO(1:50,:);
    %%
    BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
    BW_MO = BW_empty; BW_MO(indexMO) =1;
    BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
    %% prepare regressor and signal for regression
    regressor1 = zscore(VMO,[],2);
    regressor = regressor1(:,1:40000);
    Vregressor_GPU =gpuArray(regressor);

    signal1 = USSp*VSSp(:,1:40000);
    signal_GPU = gpuArray(signal1);
    %% regression
    [a1, b1, R] = CanonCor2(signal_GPU', Vregressor_GPU');
    a = gather(a1); b = gather(b1); R2 = gather(R);
    kk1 = b(:,1:50) *a(:,1:50)';
    k1_real = UMO*kk1;
    %% prepare test regressor 
    data = UselectedMO*V(:,40001:60000);    
    % subtract mean as we did before
    data = data-mean(data,2);
    regressor_test = Unew_MO' * data;   
    regressor_test_z = zscore(regressor_test,[],2);
    %%
    kernel_full2 = zeros(size(BW_empty,1)*size(BW_empty,2),size(BW_empty,1)*size(BW_empty,2));
    for j = 1:numel(indexSSp)
         kernel_temp2 = zeros(size(BW_empty,1)*size(BW_empty,2),1);
         kernel_temp2(indexMO) = k1_real(:,j);
         kernel_full2(:,indexSSp(j)) = kernel_temp2;
    end
    kernel_full2 = reshape(kernel_full2,[size(mimgT2,1),size(mimgT2,2),...
    size(mimgT2,1),size(mimgT2,2)]);
%%
    traceSSp2 = UselectedSSp*V(:,40001:60000);
    for n = 1:50
        traceSSp_predict = regressor_test_z(1:50,:)' * b(:,1:n) * a(:,1:n)';
        traceSSp_predict = traceSSp_predict';
        explained_var5(n,1) = sseExplainedCal(traceSSp2(:)',traceSSp_predict(:)');
    end
    %% color overlay a cycle
    color1 = colorcet('C06', 'N', 9);
    mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
    h1 = figure;
    axx4 = subplot(1,2,1);
    im2 = imagesc(mimgtransformed2);
    colormap(gray)
    BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
    set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');

    hold on;
    for kkk = 1:8
        TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
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
        scatter(axx4,point(kkk,2),point(kkk,1),48,...
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
    %
    ax2 = subplot(1,2,2);
    colormap(ax2,color1);
    cb = colorbar('YTickLabel', nameList,'YTick',0.06:0.13:1.05);
    cb.TickLabelInterpreter = 'none';
    axis off;
    originalSize1 = get(ax2, 'Position');
    set(ax2, 'Position', [originalSize1(1)-0.2,originalSize1(1)-0.2,originalSize1(3)*0.6, originalSize1(4)*0.3]);
    %%
    savefolder = 'AP-3/';
    print(h1, [savefolder fname '-AP'], '-dpdf', '-bestfit', '-painters');
    save([savefolder fname '-AP'],'kernel_full2','explained_var5','a','b');
    close all;
end