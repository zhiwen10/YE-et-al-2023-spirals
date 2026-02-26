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
spirals_all = [];
kk = 7;
% kk = 1:size(T,1)
clear spiralsT nframe
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');

fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','spirals_grouping',...
    [fname '_spirals_group_fftn.mat']));
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));

subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)]; 
%%
% spiral_duration = cellfun(@(x) size(x,1), archiveCell);
% indx2 = (spiral_duration>=2);
% groupedCells = archiveCell(indx2);
% filteredSpirals = cell2mat(groupedCells);   
filteredSpirals = [];
for m = 1:size(archiveCell,1)
    temp1 = archiveCell{m};
    slength = size(temp1,1);
    temp1(:,6) = slength;
    filteredSpirals = [filteredSpirals;temp1];
end
%%
U1 = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
figure;
imagesc(mimg);
%%
rate = 1;
freq = [2,8];
Fs = 35;
U2 = U(200:4:400,250:4:500,1:50);
Ur = reshape(U2, size(U2,1)*size(U2,2), size(U2,3));
meanTrace = Ur*dV(1:50,:);
meanTrace = meanTrace -mean(meanTrace ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = reshape(tracePhase,[],size(U2,1),size(U2,2));
tracePhase = permute(tracePhase,[2,3,1]);
%%
meanTrace2 = reshape(meanTrace,[],size(U2,1),size(U2,2));
%%
figure;
imagesc(mimg(200:4:400,250:4:500));
%%
figure;
% plot(1/35:1/35:400/35,squeeze(meanTrace2(1:400,35,50)));
plot(1/35:1/35:size(meanTrace2,1)/35,squeeze(meanTrace2(:,35,30)));
%%
filteredSpirals2 = filteredSpirals;
filteredSpirals2(:,1) = filteredSpirals2(:,1)-250;
filteredSpirals2(:,2) = filteredSpirals2(:,2)-200;
indx1 = (filteredSpirals2(:,1)<500-250 & filteredSpirals2(:,1)>0 ...
    & filteredSpirals2(:,2)<400-200 & filteredSpirals2(:,2)>0);
filteredSpirals2 = filteredSpirals2(indx1,:);
%%
th2 = 1:5:360; 
figure;
ax1 = subplot(1,1,1);
frame = 68666;
imagesc(squeeze(tracePhase(:,:,frame)));
colormap(ax1,colorcet('C06'));
hold on;
spiral_temp = filteredSpirals2(filteredSpirals2(:,5) == frame,:);
% scatter(filteredSpirals2(1,1)/4,filteredSpirals2(1,2)/4,24,'k','filled');
if not(isempty(spiral_temp))
    for kk = 1:size(spiral_temp,1)
        px1 = spiral_temp(kk,1)/4;
        py1 = spiral_temp(kk,2)/4;
        r = spiral_temp(kk,3)/4;
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
            color1 = 'w';
        else
            color1 = 'k';                                                      % clockwise, then color black
        end
        plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',1);          % draw the circle at max radius
        hold on;
        scatter(spiral_temp(kk,1)/4,spiral_temp(kk,2)/4,8,color1,'filled');
    end
end 
%%
cx3 = cx2(cx2<63 &cy2<51);
cy3 = cy2(cx2<63 &cy2<51);
for k = 1:length(cx3)
    tracePhaseA(k,1) = tracePhase(cy3(k),cx3(k),82);
    tracePhaseB(k,1) = tracePhase(cy3(k),cx3(k),83);
end
%%
ft = angle(exp(i*(tracePhaseB-tracePhaseA)));
%%
figure;
scatter(1:41,tracePhaseA,'k');
hold on;
scatter(1:41,tracePhaseB,'r');
hold on;
scatter(1:41,ft,'b');
%%
frt = mean(ft*35/(2*pi));