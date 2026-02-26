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
for kk = 1:size(T,1)
    clear nframe ft filteredSpirals2 meanTrace tracePhaseTemp
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)]; 
    %%
    xLim = [250,500];
    yLim = [200,400];
    xLength = xLim(2)-xLim(1);
    yLength = yLim(2)-yLim(1);
    scale1 = 4;
    freq = [2,8];
    Fs = 35;
    U2 = U(yLim(1):scale1:yLim(2),xLim(1):scale1:xLim(2),1:50);
    Ur = reshape(U2, size(U2,1)*size(U2,2), size(U2,3));
    meanTrace = Ur*dV(1:50,:);
    meanTrace = meanTrace -mean(meanTrace,2);
    meanTrace = double(meanTrace)';
    meanTrace = meanTrace(1:60000,:);                                      % only use 30 minutes data for memory
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    meanTrace = filtfilt(f1,f2,meanTrace);
    traceHilbert =hilbert(meanTrace);
    tracePhase = angle(traceHilbert);
    tracePhase = reshape(tracePhase,[],size(U2,1),size(U2,2));
    tracePhase = permute(tracePhase,[2,3,1]);
    %%
    th2 = 1:5:360; 
    px1 = round(xLength/scale1/2);
    py1 = round(yLength/scale1/2);
    r = 50/scale1;
    cx2 = round(r*cosd(th2)+px1);
    cy2 = round(r*sind(th2)+py1);
    %%
    for m = 1:numel(cy2)
        tracePhaseTemp(m,:) = tracePhase(cy2(m),cx2(m),:);
    end
    %%
    tracePhaseTemp = unwrap(tracePhaseTemp');
    instFreqAll = abs(diff(tracePhaseTemp))*35/(2*pi);
    ft = mean(instFreqAll,2);
    frt_all1{kk} = ft;
end
save('instantenous_frequency_across_sessions_dV_filt_all.mat','frt_all1');
%%
clear ratio1 N1
edges = 0:0.5:15;
for n = 1:15
    frt1 = frt_all1{n};
    counts1 = histcounts(frt1(:,1),edges);
    N1(:,n) = counts1;
    ratio1(:,n) = N1(:,n)./size(frt1,1);
end
%%
figure;
for n = 1:15
    subplot(3,5,n);
    stairs(edges(2:end),ratio1(:,n));
end
%%
frt_all1 = vertcat(frt_all1{:});
%%
ratio_mean = mean(ratio1,2);
ratio_sem = std(ratio1,[],2)./sqrt(15);
%%
h1b = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,2,1);
s1 = rectangle('Position',[2 0 6 0.1],'FaceColor',[0.5 .5 .5]);
hold on;
s2 = errorbar(edges(1:end-1),ratio_mean,ratio_sem,'k');
hold on;
[ma,mindex] = max(ratio_mean);
s3 = xline(edges(mindex),'--');
text(edges(mindex),0.08,['maxFreq = ' num2str(edges(mindex)) 'Hz']);
yticks([0:0.02:0.1]);
ylim([0,0.1]);
yticklabels({[0:0.02:0.1]});
xlabel('Instantenous frequency (Hz)');
ylabel('Proportion of spirals');
set(gca, 'box', 'off');

% subplot(1,3,2);
% frt_all2 = frt_all1(randperm(size(frt_all1,1),20000),:);
% [h,c] = scatter_kde(frt_all2(:,2), frt_all2(:,1),'filled','MarkerSize', 12);
% % Add Color bar
% cb = colorbar();
% cb.Label.String = 'Probability density estimate';
% xlim([0,15]);

subplot(1,2,2);
for n = 1:15
    clear frt
    frt = frt_all{n};
    for mm = 1:14
        clear frt_temp
        frt_temp = frt(frt(:,2)==mm,1);
        frt_mean(n,mm) = mean(frt_temp,"omitnan");
        frt_sem(n,mm) = std(frt_temp,[],1,"omitnan")./sqrt(length(frt_temp));
    end
end
frt_mean2 = mean(frt_mean,1,"omitnan");
count1 = sum(not(isnan(frt_mean)),1);
frt_std2 = std(frt_mean,[],1,"omitnan");
for m = 1:14
    frt_sem2(m,1) = frt_std2(m)./sqrt(count1(m));
end
errorbar([1:10]/35,frt_mean2(1:10),frt_sem2(1:10),'k');
xlim([0,0.3]);
ylim([2,8]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');
%%
print(h1b, 'Instantenous frequency.pdf','-dpdf', '-bestfit', '-painters');