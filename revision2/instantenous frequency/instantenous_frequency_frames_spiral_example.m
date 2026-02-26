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
pixSize = 3.45/1000/0.6*3;  % mm / pix for spiral radius
%%
kk = 7;
%%
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
% load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
load(fullfile(data_folder,'tables',[fname '_tform.mat']));                 % load atlas transformation matrix tform;
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
%%
filteredSpirals2 = cat(1,archiveCell{:});
[filteredSpirals2(:,1),filteredSpirals2(:,2)] = transformPointsForward(...
tform,filteredSpirals2(:,1),filteredSpirals2(:,2)); 
%%
frame = 58799;
spiral = filteredSpirals2(filteredSpirals2(:,5) == 58799,:);
%%
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)]; 
Ut = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
Ut2 = Ut(:,:,1:50);
Utr = reshape(Ut2,size(Ut2,1)*size(Ut2,2),size(Ut2,3));
%%
freq = [2,8];
Fs = 35;
trace = Utr*dV(1:50,frame-30:frame+30);
meanTrace = trace -mean(trace,2);
meanTrace = double(meanTrace)';             
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = reshape(tracePhase,[],size(Ut2,1),size(Ut2,2));
tracePhase = permute(tracePhase,[2,3,1]);
%%
clear cx4 cy4 cx3 cy3 tracePhase_temp
color2 = cbrewer2('seq','Greys',6);
scale3 = 5;
th2 = 1:5:360; 
th3 = 1:60:360; 
a = 63; b = 100; % height
c = 95; d = 132; % width
px1 = spiral(1,1);
py1 = spiral(1,2);
r = spiral(1,3)*pixSize/(0.01);
%%
clear ft ft2 tracePhaseA tracePhaseB
th2 = 1:5:360; 
cx2 = []; cy2 = [];
for r1 = 16:16:r
    cx2 = [cx2,round(r1*cosd(th2)+px1)];
    cy2 = [cy2,round(r1*sind(th2)+py1)];
end
for k = 1:length(cx2)
    tracePhaseA(k,1) = squeeze(tracePhase(cy2(k),cx2(k),30));
    tracePhaseB(k,1) = squeeze(tracePhase(cy2(k),cx2(k),31));
end
ft = angle(mean(exp(i*(tracePhaseB-tracePhaseA))));
ft = abs(ft);
ft2 = ft*35/(2*pi);
%%
BW2a = imresize(BW2,[1320,1140]);
ra  = [30:16:r];
h1a = figure('Renderer', 'painters', 'Position', [100 100 400 300]);
for k = 1:2
    ax1 = subplot(2,2,2*k-1);
    im_phase = imagesc(tracePhase(:,:,29+k));
    set(im_phase, 'AlphaData', BW2a, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([-pi,pi]);
    colormap(ax1,colorcet('C06'))
    axis image; axis off;
    count = 1;
    for r1 = ra
        cx3 = r1*cosd(th2)+px1;
        cy3 = r1*sind(th2)+py1;
        cx4 = r1*cosd(th3)+px1;
        cy4 = r1*sind(th3)+py1;
        hold on;
        if spiral(1,4) == 1                                               % counterclockwise, then color white
            color1 = 'w';
        else
            color1 = 'k';                                                      % clockwise, then color black
        end
        plot([cx3 cx3(1)],[cy3 cy3(1)],'color',color2(count+2,:),'LineWidth',1);          % draw the circle at max radius
        hold on;
        scatter(px1,py1,8,'filled','MarkerFaceColor',color2(count+2,:));
        hold on;
        scatter(cx4,cy4,8,'filled','MarkerFaceColor',color2(count+2,:));
        count = count+1;
    end
    ylim([a,b]*8);
    xlim([c,d]*8);
end
%
ax2 = subplot(3,2,2);
r1 = ra(end);
cx4 = round(r1*cosd(th3)+px1);
cy4 = round(r1*sind(th3)+py1);
for m = 1:numel(cx4)
    for k = 1:2
        tracePhase_temp(m,k) = tracePhase(cy4(m),cx4(m),29+k);
    end
end
tracePhase_temp = unwrap(tracePhase_temp,[],1);
plot(1:numel(cx4),tracePhase_temp(:,1));
hold on;
plot(1:numel(cx4),tracePhase_temp(:,2));
phase_diff = abs(tracePhase_temp(:,2)-tracePhase_temp(:,1));    
xlim([1,6]);
xticks([1:6]);
box off;
count = 1;
for r1 = ra(end)
    cx4 = round(r1*cosd(th3)+px1);
    cy4 = round(r1*sind(th3)+py1);
    for m = 1:numel(cx4)
        for k = 1:2
            tracePhase_temp(m,k) = tracePhase(cy4(m),cx4(m),29+k);
        end
    end
    tracePhase_temp = unwrap(tracePhase_temp,[],1);
    phase_diff = abs(tracePhase_temp(:,2)-tracePhase_temp(:,1));  
%     ax4 = subplot(3,2,4);
%     plot(1:numel(cx4),phase_diff,'color',color2(6,:));
%     hold on;
%     ylim([0.4,1.2]);
%     box off;
    ax6 = subplot(3,2,4);
    freq1 = phase_diff*35/(2*pi);
    plot(1:numel(cx4),freq1,'color',color2(6,:));
    hold on;
    xlim([1,6]);
    xticks([1:6]);
    ylim([2,8]);
    box off;
    count = count+1;
end
count = 1;
for r1 = ra
    cx4 = round(r1*cosd(th3)+px1);
    cy4 = round(r1*sind(th3)+py1);
    for m = 1:numel(cx4)
        for k = 1:2
            tracePhase_temp(m,k) = tracePhase(cy4(m),cx4(m),29+k);
        end
    end
    tracePhase_temp = unwrap(tracePhase_temp,[],1);
    phase_diff = abs(tracePhase_temp(:,2)-tracePhase_temp(:,1));  
    ax6 = subplot(3,2,6);
    freq1 = phase_diff*35/(2*pi);
    plot(1:numel(cx4),freq1,'color',color2(count+2,:));
    hold on;
    xlim([1,6]);
    xticks([1:6]);
    ylim([2,8]);
    box off;
    count = count+1;
end
ax6 = subplot(3,2,6);
yline(ft2,'r--');
print(h1a, 'Instantenous frequency_example.pdf','-dpdf', '-bestfit', '-painters');  