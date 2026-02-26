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
for kk = 1:size(T,1)
    clear nframe frt filteredSpirals2 meanTrace
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
%     filteredSpirals = [];
%     for m = 1:size(archiveCell,1)
%         temp1 = archiveCell{m};
%         slength = size(temp1,1);
%         temp1(:,6) = slength;
%         filteredSpirals = [filteredSpirals;temp1];
%     end
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    archiveCell = archiveCell(indx2);
    filteredSpirals2 = cell2mat(archiveCell);  
    % filteredSpirals2 = filteredSpirals(filteredSpirals(:,5)<60000,:);
    % filteredSpirals2 = filteredSpirals(filteredSpirals(:,6)>8,:);
    [filteredSpirals2(:,1),filteredSpirals2(:,2)] = transformPointsForward(...
    tform,filteredSpirals2(:,1),filteredSpirals2(:,2));  
    filteredSpirals2(:,1:2) = round(filteredSpirals2(:,1:2)); 
    indx1 = (filteredSpirals2(:,1)<950 & filteredSpirals2(:,1)>850 ...
        & filteredSpirals2(:,2)<700 & filteredSpirals2(:,2)>500 & filteredSpirals2(:,3)>=80);
    filteredSpirals2 = filteredSpirals2(indx1,:);
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)]; 
    %%
    scale3 = 5;
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    h1ac = figure('Renderer', 'painters', 'Position', [100 100 700 700]);
    ax1 = subplot(1,1,1);
    im1= imshow(mimgt);
    set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    lineColor = 'k'; lineColor1 = 'w';
    hemi = [];
    overlayOutlines(coords,1);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis on;
    %%
    Ut = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    %%
    freq = [2,8];
    Fs = 35;
    for m = 1:7
        trace = squeeze(Ut(pixel(m,1),pixel(m,2),1:50))'*dV(1:50,:);
        meanTrace = trace -mean(trace);
        meanTrace = double(meanTrace)';             
        [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
        meanTrace = filtfilt(f1,f2,meanTrace);
        traceHilbert =hilbert(meanTrace);
        tracePhase = angle(traceHilbert);
        %%
        for j = 1:size(filteredSpirals2,1)
            phase_all(j,m) = tracePhase(filteredSpirals2(j,5));
            phase_diff = tracePhase(filteredSpirals2(j,5)+1)-tracePhase(filteredSpirals2(j,5));
            freq_all(j,m) = angle(exp(i*(phase_diff)))*35./(2*pi);
        end
    end
    %%
    phase_bin = -pi:pi/6:pi;
    figure;
    for  m = 1:7
        subplot(7,1,m);
        phase_temp = phase_all(:,m);
        freq_temp = freq_all(:,m);    
        %%
        for bini = 1:numel(phase_bin)-1
            index = find(phase_temp>phase_bin(bini) & phase_temp<phase_bin(bini+1));
            freq_mean(bini,1) = mean(freq_temp(index));
            % freq_sem(bini,1) = std(freq_temp(index))./sqrt(sum(index));
            freq_sem(bini,1) = std(freq_temp(index));
        end
        scatter(phase_temp,freq_temp,4,'k');   
        hold on;
        errorbar(phase_bin(1:end-1),freq_mean,freq_sem,'r')
        ylim([0,8]);
    end
end
%%
save('instantenous_frequency_across_sessions_dV_filt.mat','frt_all');
%%
edges = 0:0.5:15;
for n = 1:15
    frt1 = frt_all{n};
    counts1 = histcounts(frt1(:,1),edges);
    N1(:,n) = counts1;
    ratio1(:,n) = N1(:,n)./size(frt1,1);
end
%%
figure;
for n = 1:15
    subplot(3,5,n);
    stairs(edges(2:end),log10(ratio1(:,n)));
end
%%
frt_all1 = vertcat(frt_all{:});
%%
ratio_mean = mean(ratio1,2);
ratio_sem = std(ratio1,[],2)./sqrt(15);
%%
h1b = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,2,1);
s1 = rectangle('Position',[2 0 6 0.25],'FaceColor',[0.5 .5 .5]);
hold on;
s2 = errorbar(edges(1:end-1),ratio_mean,ratio_sem,'k');
hold on;
[ma,mindex] = max(ratio_mean);
s3 = xline(edges(mindex),'--');
text(edges(mindex),0.08,['maxFreq = ' num2str(edges(mindex)) 'Hz']);
yticks([0:0.05:0.25]);
ylim([0,0.25]);
yticklabels({[0:0.05:0.25]});
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
ylim([4,5]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');
%%
print(h1b, 'Instantenous frequency_filtered.pdf','-dpdf', '-bestfit', '-painters');