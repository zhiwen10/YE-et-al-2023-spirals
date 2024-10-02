function getAlphaSprialMap(T,freq,data_folder,save_folder)
freq_name = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz']; 
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
%% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
xx = 1:1140; yy = 1:1320;
[xxq,yyq] = meshgrid(xx,yy);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
%%
hist_bin = 40;
pixSize = 0.01;                                                            % pixel resolution mm/pix
pixArea = pixSize^2;                                                       % 2d-pixel resolution mm^2/pix^2
%%
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
pixel(8,:) = [290,700]; % MOs
pixel = round(pixel/params.downscale);
color1 = cbrewer2('qual','Set3',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
%%
for kk = 1:size(T,1)
    clear a1 b1 a2 b2 filteredSpirals_alpha filteredSpirals_nonalpha...
    unique_spirals_alpha unique_spirals_nonalpha...
    unique_spirals_unit_alpha unique_spirals_unit_nonalpha...
    spiralsT filteredSpirals rawTrace traceAmp1 tracePhase1...
    tracefilt1 traceAmp traceraw tracefilt amp traceAlpha1 traceNonalpha1
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals\svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    dV = dV(1:50,:);
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));
    load(fullfile(data_folder,'spirals\spirals_freq\spirals_fftn\',freq_name,...
        [fname '_spirals_group_fftn.mat']));
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));    
    filteredSpirals(:,1:2) = round(spiralsT);
    %%
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % find spirals within the brain boundry
    filteredSpirals = filteredSpirals(lia,:);
    %%
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,1:50);
    mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
    %%
    for i = 1:8
        rawTrace(i,:) = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'*V(1:50,:);
        rawTrace(i,:) = rawTrace(i,:)./mimgtransformed(pixel(i,1),pixel(i,2));
    end
    %%
    Fs = 35;
    rawTrace = double(rawTrace);
    rawTrace = rawTrace -mean(rawTrace ,2); 
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    traceFilt = filtfilt(f1,f2,rawTrace');
    traceFilt = traceFilt';
    traceHilbert =hilbert(traceFilt);
    tracePhase = angle(traceHilbert);
    traceAmp = abs(traceHilbert);
    %%
    epochLength = 70;
    [Tall,Ta,Tna,thred1,thred2,traceAlpha,traceNonalpha] = getAmpEpochs3(traceAmp,epochLength,rawTrace);
    %%
    traceAlpha = traceAlpha(3:7,:,:);
    traceNonalpha = traceNonalpha(3:7,:,:);
    traceAlpha1 = reshape(traceAlpha,size(traceAlpha,1)*size(traceAlpha,2),[]);
    traceNonalpha1 = reshape(traceNonalpha,size(traceNonalpha,1)*size(traceNonalpha,2),[]);
    [freq1,psdx_alpha,psdx_mean] = fft_spectrum(traceAlpha1);
    [freq1,psdx_nonalpha,non_psdx_mean] = fft_spectrum(traceNonalpha1);
%     psdx_alpha = reshape(psdx_alpha,size(traceNonalpha,1),size(traceNonalpha,2),[]);
%     psdx_nonalpha = reshape(psdx_nonalpha,size(traceNonalpha,1),size(traceNonalpha,2),[]);
    psdx_alpha_mean = mean(psdx_alpha,2);
    psdx_nonalpha_mean = mean(psdx_nonalpha,2);
    %%
    h1 = figure('Renderer', 'painters', 'Position', [100 100 250 400]);
    histogram(Tall.traceAmp);
    hold on;
    xline(thred1,'r');
    xline(thred2,'r');
    print(h1, fullfile(save_folder,[fname '_amp_historgam']),...
    '-dpdf', '-bestfit', '-painters');
    %%
    h2 = figure('Renderer', 'painters', 'Position', [100 100 900 500]);
    for i = 3:7
        plot(t,traceFilt(i,:)/5+i*0.01,'color',color1(i,:));
        hold on;
    end
    stairs(t(Tall.firstFrame),Tall.traceAmp,'r');
    hold on;
    for k = 1:size(Ta.firstFrame,1)
        plot([t(Ta.firstFrame(k)),t(Ta.lastFrame(k))]',...
            [Ta.traceAmp(k),Ta.traceAmp(k)]','b','lineWidth',7);
    end
    hold on;
    for k = 1:size(Tna.firstFrame,1)
        plot([t(Tna.firstFrame(k)),t(Tna.lastFrame(k))]',...
            [Tna.traceAmp(k),Tna.traceAmp(k)]','g','lineWidth',7);
    end
    print(h2, fullfile(save_folder,[fname '_traces']),...
    '-dpdf', '-bestfit', '-painters');
    %%
    h3 = figure('Renderer', 'painters', 'Position', [100 100 250 400]);
    freq_value = [0.5,2,4,6,8];
    log_freq_value = log(freq_value);
    freqN = 17;
    plot(log(freq1(2:freqN)),log(psdx_alpha_mean(2:freqN)),'r');
    hold on;
    plot(log(freq1(2:freqN)),log(psdx_nonalpha_mean(2:freqN)),'k');
    xlim([log_freq_value(1),log_freq_value(end)]);
    xticks(log_freq_value);
    xticklabels({'0.5','2','4','6','8'});
    xlabel('log Frequency');
    ylabel('log Power (df/f^2)');
    print(h3, fullfile(save_folder,[fname '_fft_alpha']),...
    '-dpdf', '-bestfit', '-painters');
    %% get frames for alpha and non-alpha
    alphaFrames = [];
    for i = 1:size(Ta,1)
        temp = [Ta.firstFrame(i):Ta.lastFrame(i)]';
        alphaFrames = [alphaFrames;temp];
    end
    nonalphaFrames = [];
    for i = 1:size(Tna,1)
        temp = [Tna.firstFrame(i):Tna.lastFrame(i)]';
        nonalphaFrames = [nonalphaFrames;temp];
    end
    %% get sprials in alpha and non-alpha frames    
    [a1,b1] = ismember(filteredSpirals(:,5),alphaFrames);
    filteredSpirals_alpha = filteredSpirals(a1,:);
    [a2,b2] = ismember(filteredSpirals(:,5),nonalphaFrames);
    filteredSpirals_nonalpha = filteredSpirals(a2,:);
    %% get spiral density map for alpha and non-alpha frames
    [unique_spirals_alpha,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(filteredSpirals_alpha,hist_bin);
    [unique_spirals_nonalpha,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(filteredSpirals_nonalpha,hist_bin);
    %% desnity line for spirals in alpha frames
    nframe = numel(alphaFrames);
    unique_spirals_unit_alpha = unique_spirals_alpha(:,3)/...
        (hist_bin*hist_bin*pixArea);                                           % spiral counts/mm^2
    unique_spirals_unit_alpha = unique_spirals_unit_alpha./nframe*35;                      % spirals/(mm^2*s)
    % interp histgram counts 
    F = scatteredInterpolant(unique_spirals_alpha(:,1),...
        unique_spirals_alpha(:,2),...
        unique_spirals_unit_alpha);
    vq1 = F(xxq,yyq);
    for i = 1:size(points,1)
        count_sample_alpha(i,1) = vq1(points(i,2),points(i,1));
    end
    %% desnity line for spirals in non-alpha frames
    nframe_na = numel(nonalphaFrames);
    unique_spirals_unit_nonalpha = unique_spirals_nonalpha(:,3)/...
        (hist_bin*hist_bin*pixArea);                                           % spiral counts/mm^2
    unique_spirals_unit_nonalpha = unique_spirals_unit_nonalpha./nframe_na*35;                      % spirals/(mm^2*s)
    % interp histgram counts 
    Fna = scatteredInterpolant(unique_spirals_nonalpha(:,1),...
        unique_spirals_nonalpha(:,2),...
        unique_spirals_unit_nonalpha);
    vq1_na = Fna(xxq,yyq);
    for i = 1:size(points,1)
        count_sample_nonalpha(i,1) = vq1_na(points(i,2),points(i,1));
    end
    %%
    h4 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
    % plot spirals density map
    ax1 = subplot(1,3,1);
    scatter(ax1,unique_spirals_alpha(:,1),unique_spirals_alpha(:,2),...
        1,unique_spirals_unit_alpha,'filled');                                   % plot scatter plot with color indicating spiral density 
    max_c = 3;
    caxis([0,max_c]);
    colormap(ax1,hot);
    cb0 = colorbar(ax1); 
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    axis off; axis image;
    xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
    
    cb0.Ticks = [0,max_c];
    cb0.TickLabels = {'0',num2str(max_c)};
    
    ax2 = subplot(1,3,2);
    scatter(ax2,unique_spirals_nonalpha(:,1),unique_spirals_nonalpha(:,2),...
        1,unique_spirals_unit_nonalpha,'filled');                                   % plot scatter plot with color indicating spiral density 
    caxis([0,max_c]);
    colormap(ax2,hot);
    cb1 = colorbar(ax2); 
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    axis off; axis image;
    xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
    
    cb1.Ticks = [0,max_c];
    cb1.TickLabels = {'0',num2str(max_c)};
    
    subplot(1,3,3)
    points1 = points-points(1,:);
    points_line = vecnorm(points1,2,2);
    plot(points_line(:,1),count_sample_alpha,'r');
    hold on;
    plot(points_line(:,1),count_sample_nonalpha,'k');
    xlim([0,800]);
    ylim([0,5]);
    yticks([0,0.5,1,1.5,2]);
    yticklabels({'0' '0.5' '1' '1.5' '2'});
    xticks([0 200 400 600 800]);
    xticklabels({'0','2','4','6','8'});
    ylabel('sprials/(mm^2*s)'); 
    legend({'Alpha','Non-alpha'});
    print(h4, fullfile(save_folder,[fname '_density_map_alpha']),...
    '-dpdf', '-bestfit', '-painters');
    %%
    close all;
    %%
    save(fullfile(save_folder,[fname '_histogram.mat']),...
        'unique_spirals_alpha','unique_spirals_nonalpha',...
        'alphaFrames','nonalphaFrames',...
        'count_sample_alpha','count_sample_nonalpha',...
        'unique_spirals_unit_alpha','unique_spirals_unit_nonalpha',...
        'filteredSpirals_alpha','filteredSpirals_nonalpha');
end
