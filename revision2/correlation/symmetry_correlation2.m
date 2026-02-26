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
scale = 8;
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
trace_raw_all = nan(60000,2,15);
trace_raw_hemi = nan(60000,2,15);
trace_raw_ap = nan(60000,2,15);
kkk = [1,8];
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
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    %%
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    %%
%     V = [zeros(size(V,1),1) diff(V,[],2)];
%     V = double(V);
%     Fs = 35;
%     [f1,f2] = butter(2, [2,8]/(Fs/2), 'bandpass');
%     V = filtfilt(f1,f2,V');
%     V = V';
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    %% hemi
    load(fullfile(data_folder,'spirals_mirror','regression_kernels',...
    'kernelMaps_allSessions_hemi.mat'));
    for m = 1:numel(kkk)
        TheColorImage = squeeze(colorIntensityAll(:,:,kkk(m),kk));  
        maxI = max(TheColorImage(:));
        [row,col] = find(TheColorImage==maxI);
        trace_raw1 = squeeze(Utransformed(point(kkk(m),1)*8,point(kkk(m),2)*8,1:50))'* V(1:50,1:60000);
        trace_raw_mean1 = mimgtransformed(point(kkk(m),1)*8,point(kkk(m),2)*8); %sensory
        trace_raw_all1(:,m,kk) = trace_raw1./trace_raw_mean1;
        trace_raw2 = squeeze(Utransformed(row*scale,col*scale,1:50))'* V(1:50,1:60000);
        trace_raw_mean2 = mimgtransformed(row*8,col*8); %left hemisphere
        trace_raw_hemi(:,m,kk) = trace_raw2./trace_raw_mean2;
    end
    %% ap
    load(fullfile(data_folder,'spirals_mirror','regression_kernels',...
    'kernelMaps_allSessions_AP.mat'));
    for m = 1:numel(kkk)
        TheColorImage = squeeze(colorIntensityAll(:,:,kkk(m),kk));  
        maxI = max(TheColorImage(:));
        [row,col] = find(TheColorImage==maxI);
        trace_raw1 = squeeze(Utransformed(point(kkk(m),1)*8,point(kkk(m),2)*8,1:50))'* V(1:50,1:60000);
        trace_raw_mean1 = mimgtransformed(point(kkk(m),1)*8,point(kkk(m),2)*8); %sensory
        trace_raw_all2(:,m,kk) = trace_raw1./trace_raw_mean1;
        trace_raw2 = squeeze(Utransformed(row*scale,col*scale,1:50))'* V(1:50,1:60000);
        trace_raw_mean2 = mimgtransformed(row*8,col*8); %left hemisphere
        trace_raw_ap(:,m,kk) = trace_raw2./trace_raw_mean2;
    end
end
%%
tsteps = 70;
lagN = tsteps*2+1;
% lag * pixel * mouse * corr_pair
c1 = nan(lagN,2,15,3);
for m = 1:2
    for i = 1:15
        [c1(:,m,i,1),lags] = xcorr(trace_raw_all1(:,m,i),trace_raw_hemi(:,m,i),tsteps,'normalized');
    end
end

for m = 1:2
    for i = 1:15
        [c1(:,m,i,2),lags] = xcorr(trace_raw_all1(:,m,i),trace_raw_ap(:,m,i),tsteps,'normalized');
    end
end


for m = 1:2
    for i = 1:15
        [c1(:,m,i,3),lags] = xcorr(trace_raw_all1(:,m,i),trace_raw_all1(:,m,i),tsteps,'normalized');
    end
end
%%
c1_ssp = squeeze(c1(:,1,:,:));
c1_mean = squeeze(mean(c1_ssp,2,'omitnan'));
c1_sem = squeeze(std(c1_ssp,[],2,'omitnan')/sqrt(15));
%%
pair_names = {'HEMI','AP','AUTO'};
color1 = {'r','g','b'};
hs = figure('Renderer', 'painters', 'Position', [100 100 900 250]);
tlim = 2;
% ticks = [-0.5:0.25:0.5];
ticks = [-2:1:2];
for m  = 1:3
    subplot(1,3,m);
    for i = 1:15
        plot(lags/35,c1_ssp(:,i,m),'k');
        hold on;
    end
    shadedErrorBar(lags/35, c1_mean(:,m), c1_sem(:,m), 'lineprops', ['-' color1{m}]);
    yline(0,'k--');
    hold on;
    xline(0,'k--');
    xlim([-tlim,tlim]);
    ylim([-0.5,1]);
    % xticks([-2:1:2]);
    % xticklabels([-2:1:2]);
    xticks(ticks);
    xticklabels(ticks);
    title(pair_names{m});
    xlabel('Time lag (s)');
    ylabel('Correlation coefficient');
end
%%
print(hs, 'correlation_broad_spectrum.pdf','-dpdf', '-bestfit', '-painters');
% print(hs, 'correlation_2_8Hz.pdf','-dpdf', '-bestfit', '-painters');