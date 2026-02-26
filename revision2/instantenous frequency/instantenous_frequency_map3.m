githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
%%
freq = [2,8];
Fs = 35;
freq_all1 = [];
for kk = 1:size(T,1)
    clear nframe frt filteredSpirals2 meanTrace
    freq_all = [];
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    archiveCell = archiveCell(indx2);
    filteredSpirals2 = cell2mat(archiveCell);  
    filteredSpirals2(filteredSpirals2(:,4)== 0,4) = -1;
    filteredSpirals2 = filteredSpirals2(filteredSpirals2(:,5)<60000,:);
    [filteredSpirals2(:,1),filteredSpirals2(:,2)] = transformPointsForward(...
    tform,filteredSpirals2(:,1),filteredSpirals2(:,2));  
    filteredSpirals2(:,1:2) = round(filteredSpirals2(:,1:2)); 
    indx1 = (filteredSpirals2(:,1)<950 & filteredSpirals2(:,1)>850 ...
        & filteredSpirals2(:,2)<700 & filteredSpirals2(:,2)>500 ...
        & filteredSpirals2(:,3)>=80 & filteredSpirals2(:,4)== -1);
%     indx1 = (filteredSpirals2(:,1)<950 & filteredSpirals2(:,1)>850 ...
%         & filteredSpirals2(:,2)<700 & filteredSpirals2(:,2)>500 ...
%         & filteredSpirals2(:,3)>=80);
    filteredSpirals2 = filteredSpirals2(indx1,:);
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)]; 
    %%
%     scale3 = 5;
%     mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%     h1ac = figure('Renderer', 'painters', 'Position', [100 100 700 700]);
%     ax1 = subplot(1,1,1);
%     im1= imshow(mimgt);
%     set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
%     hold on;
%     lineColor = 'k'; lineColor1 = 'w';
%     hemi = [];
%     overlayOutlines(coords,1);
%     hold on;
%     plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
%     axis on;
    %%
    Ut = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Ut2 = Ut(1:16:end,1:16:end,1:50);
    Utr = reshape(Ut2,size(Ut2,1)*size(Ut2,2),size(Ut2,3));
    %%
    trace = Utr*dV(1:50,1:60000);
    meanTrace = trace -mean(trace,2);
    meanTrace = double(meanTrace)';             
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    meanTrace = filtfilt(f1,f2,meanTrace);
    traceHilbert =hilbert(meanTrace);
    tracePhase = angle(traceHilbert);
    tracePhase = reshape(tracePhase,[],size(Ut2,1),size(Ut2,2));
    tracePhase = permute(tracePhase,[2,3,1]);
    for m = 1:size(filteredSpirals2,1)
            tracePhase1 = squeeze(tracePhase(:,:,filteredSpirals2(m,5)));
            tracePhase2 = squeeze(tracePhase(:,:,filteredSpirals2(m,5)+1));
            phase_diff = tracePhase2-tracePhase1;
            freq_all(:,:,m) = angle(exp(i*(phase_diff)))*35./(2*pi);
    end
    freq_all1{kk} = freq_all;
end
%%
freq_all2 = cellfun(@(x) mean(x,3), freq_all1,'UniformOutput',false);
freq_all3 = cat(3,freq_all2{:});
freq_mean = mean(freq_all3,3);
pixel2 = round(pixel./16);
for kk = 1:size(pixel2,1)
    freq1(kk,:) = squeeze(freq_all3(pixel2(kk,1),pixel2(kk,2),:));
end
freq1_mean = mean(freq1,2);
freq1_sem = std(freq1,[],2)./sqrt(15);
freq1_mean_all = mean(freq1(:));
freq_mean2 = imresize(freq_mean,[1320,1140]);
%%
color1 = cbrewer2('div','RdBu',20);
color1 = flipud(color1);
color2 = cbrewer2('seq','Greys',9);
h1ac = figure('Renderer', 'painters', 'Position', [100 100 700 700]);
ax1 = subplot(1,1,1);
im1= imagesc(freq_mean2);
colormap(ax1,color1);
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
overlayOutlines(coords,1);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(3:end,:),'filled');
caxis([3,5]);
colorbar;
axis image; axis off;
roi = drawpolygon;
bw = createMask(roi);
%%
save('insta_freq_roi_mask.mat','bw');
%%
color1 = cbrewer2('div','RdBu',20);
color1 = flipud(color1);
color2 = cbrewer2('seq','Greys',11);
load('insta_freq_roi_mask.mat');
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 300]);
ax1 = subplot(1,2,1);
im1= imagesc(freq_mean2);
colormap(ax1,color1);
BW2 = BW&bw;
set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5;
overlayOutlines(coords,1);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(5:end,:),'filled');
axis image; axis off;
caxis([3,5]);
colorbar;
ax2 = subplot(1,2,2);
for kk = 1:7
    scatter(kk,freq1(kk,:),6,'k');
    hold on;
end
yline(freq1_mean_all,'k--');
errorbar(1:7,freq1_mean,freq1_sem,'r');
area_names = {'VISp','RSP','SSp-ul','SSp-ll','SSp-m','SSp-n','SSp-bfd'};
xticks([1:7]);
xticklabels(area_names);
ylabel('Instananeous frequency (Hz)');
%%
print(h1, 'instananeous frequency map','-dpdf', '-bestfit', '-painters');
%% N-way anova
subj = [1:15];
subject = repmat(subj,[7,1]);
area = [1:7]';
areas = repmat(area,[1,15]);
subject = categorical(subject(:));
areas = categorical(areas(:));
freq2 = freq1(:);
[pp_freq,tbl,stats] = anovan(freq2,{subject,areas},...
    'model','interaction','random',1,'varnames',{'subj','area'});
%% multiple comparison
pairs3 = nchoosek(1:7, 2);  % Choose 2 elements from 1:7
for k = 1:size(pairs3,1)
    [hh1(k,1),pp1(k,1)] = ttest(freq1(pairs3(k,1),:),freq1(pairs3(k,2),:));
end
[corrected_p, hh2]=bonf_holm(pp1,0.05);
%%
save('instantenous_frequency_across_sessions_dV_filt.mat','frt_all');
