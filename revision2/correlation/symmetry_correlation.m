%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
point(8,:) = [115 105]; %% VISp
point(7,:) = [97 79]; %% RSP
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% mask and Kernel regression map for right and left
spath = string(st.structure_id_path);
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
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
% mask and Kernel regression map for left
hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);

BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
BW_right = BW_empty; BW_right(indexright) =1;
%%
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'tables',[fname '_tform.mat']));               % load atlas transformation matrix tform;
    %%
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    %%
    load(fullfile(data_folder,'spirals_mirror','regression_kernels',...
    'kernelMaps_allSessions_hemi.mat'));
    %% color overlay a cycle
%     color1 = colorcet('C06', 'N', 9);
%     figure;
%     scale = 8;
%     mimgtransformed2 = squeeze(mimgtransformed2All(:,:,kk));
%     axx4 = subplot(1,1,1);
%     im2 = imagesc(mimgtransformed2);
%     colormap(gray)
%     BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
%     set(im2, 'AlphaData', BW_left, 'AlphaDataMapping', 'scaled');
% 
%     hold on;
%     for kkk = 1:8
%         colorIntensity = squeeze(colorIntensityAll(:,:,kkk,kk));
%         imageSize = size(colorIntensity); 
%         thisColor = color1(kkk,:);
%         thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
%             thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
%         hold on
%         hb2(kkk) = imshow(thisColorImage); 
%         set(hb2(kkk),'AlphaData',colorIntensity);
%         axis image; 
%     end
%     hold on;
%     scale3 = 5/8;
%     lineColor = 'k'; lineColor1 = 'w';
%     hemi = 'left';
%     hold on;
%     plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
%     plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
%     plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
%     plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
%     set(gca,'Ydir','reverse')
%     axis image; axis off
%     kkk = [1,8];
%     for m = 1:numel(kkk)
%         TheColorImage = squeeze(colorIntensityAll(:,:,kkk(m),kk));
%         % overlay kernel with colormap    
%         maxI = max(TheColorImage(:));
%         [row,col] = find(TheColorImage==maxI);
%         hold on;
%         scatter(col,row,12,'MarkerFaceColor','k','MarkerEdgeColor','None');
%     end
    %%
    for m = 1:numel(kkk)
        TheColorImage = squeeze(colorIntensityAll(:,:,kkk(m),kk));  
        maxI = max(TheColorImage(:));
        [row,col] = find(TheColorImage==maxI);
        trace_raw1 = squeeze(Utransformed(point(kkk(m),1)*8,point(kkk(m),2)*8,1:50))'* V(1:50,1:60000);
        trace_raw_mean1 = mimgtransformed(point(kkk(m),1)*8,point(kkk(m),2)*8); %sensory
        trace_raw_all(:,m) = trace_raw1./trace_raw_mean1;
        trace_raw2 = squeeze(Utransformed(row*scale,col*scale,1:50))'* V(1:50,1:60000);
        trace_raw_mean2 = mimgtransformed(row*8,col*8); %left hemisphere
        trace_raw_hemi(:,m) = trace_raw2./trace_raw_mean2;
    end
%     %%
%     figure;
%     plot(t(1:60000),trace_raw_all(:,1),'k');
%     hold on;
%     plot(t(1:60000),trace_raw_hemi(:,1),'r');
%     
%     hold on;
%     plot(t(1:60000),trace_raw_all(:,2)+0.1,'k');
%     hold on;
%     plot(t(1:60000),trace_raw_hemi(:,2)+0.1,'r');
%     %%
%     figure;
%     subplot(1,2,1);
%     [c,lags] = xcorr(trace_raw_all(:,1),trace_raw_hemi(:,1),70,'normalized');
%     plot(lags/35,c);
%     ylim([0,1]);
%     subplot(1,2,2);
%     [c,lags] = xcorr(trace_raw_all(:,2),trace_raw_hemi(:,2),70,'normalized');
%     plot(lags/35,c);
%     ylim([0,1]);
end