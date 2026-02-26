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
figure;
imagesc(mimg);
%%
rate = 1;
freq = [2,8];
Fs = 35;
U2 = U(200:4:400,250:4:500,1:50)./mimg(200:4:400,250:4:500);
Ur = reshape(U2, size(U2,1)*size(U2,2), size(U2,3));
%% raw phase v
meanTraceV = Ur*V(1:50,:);
meanTraceV = meanTraceV -mean(meanTraceV ,2);
meanTraceV = double(meanTraceV)';
% traceHilbert =hilbert(meanTraceV);
% tracePhase = angle(traceHilbert);
% tracePhase = reshape(tracePhase,[],size(U2,1),size(U2,2));
% tracePhase = permute(tracePhase,[2,3,1]);
meanTrace_V = reshape(meanTraceV,[],size(U2,1),size(U2,2));
%% raw phase dv
meanTrace = Ur*dV(1:50,:);
meanTrace = meanTrace -mean(meanTrace ,2);
meanTrace = double(meanTrace)';
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = reshape(tracePhase,[],size(U2,1),size(U2,2));
tracePhase = permute(tracePhase,[2,3,1]);
meanTrace_raw = reshape(meanTrace,[],size(U2,1),size(U2,2));
%%
[f1,f2] = butter(3,17/(Fs/2), 'low');
meanTrace_low = filtfilt(f1,f2,meanTrace);
meanTrace_low2 = reshape(meanTrace_low,[],size(U2,1),size(U2,2));
clear meanTrace_low;
%% filtered phase
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
meanTrace1 = filtfilt(f1,f2,meanTrace);
traceHilbert1 =hilbert(meanTrace1);
tracePhase1 = angle(traceHilbert1);
clear traceHilbert1
tracePhase1 = reshape(tracePhase1,[],size(U2,1),size(U2,2));
tracePhase1 = permute(tracePhase1,[2,3,1]);
meanTrace_filt = reshape(meanTrace1,[],size(U2,1),size(U2,2));
%%
figure;
imagesc(mimg(200:4:400,250:4:500));
%%
pixel = [200+35*4,250+30*4];
% epoch = [1955,1970];
epoch = [1681,1686];
frameStart = find(t>epoch(1),1,'first'); frameEnd = find(t>epoch(2),1,'first');
sample_epoch = frameStart:frameEnd;
% sample_epoch = epoch(1)*35:epoch(2)*35;
%%
tracePhaseC = squeeze(tracePhase(35,30,sample_epoch));
tracePhaseC = unwrap(tracePhaseC);
instFreq = abs(diff(tracePhaseC))*35/(2*pi);
instFreqSmooth = smoothdata(instFreq,"gaussian",17);
%% trace example
figure;
subplot(4,1,1);
plot(sample_epoch(1)/35:1/35:sample_epoch(end)/35,squeeze(meanTrace_V(sample_epoch,35,30)));
subplot(4,1,2);
plot(sample_epoch(1)/35:1/35:sample_epoch(end)/35,squeeze(meanTrace_raw(sample_epoch,35,30)));
subplot(4,1,3);
plot(sample_epoch(1)/35:1/35:sample_epoch(end)/35,squeeze(meanTrace_filt(sample_epoch,35,30)));
subplot(4,1,4);
plot(sample_epoch(1)/35:1/35:sample_epoch(end-1)/35,instFreq,'r');
hold on;
plot(sample_epoch(1)/35:1/35:sample_epoch(end-1)/35,instFreqSmooth,'k');
%%
meanTrace2 = Ur*dV(1:50,:);
meanTrace2 = meanTrace2 -mean(meanTrace2 ,2);
meanTrace2 = double(meanTrace2)';
meanTrace_raw2 = reshape(meanTrace2,[],size(U2,1),size(U2,2));
meanTrace_raw3 = squeeze(meanTrace_raw2(:,35,30));
%%
meanTrace_raw1 = squeeze(meanTrace_raw(:,15,30));
meanTrace_filt1 = squeeze(meanTrace_filt(:,15,30));
Fs = 35;
epochSize = 20*Fs;
batchN = floor(size(meanTrace_raw1,1)/epochSize);
clear dataMatrix3
figure;
subplot(1,1,1);
% dataMatrix = zeros(batchN, epochSize);
% for i = 1:batchN
%     dataMatrix(i,:) =  meanTrace_raw1(1+epochSize*(i-1):epochSize*i);
% end   
% alpha_trace_mean = bsxfun(@minus, dataMatrix, mean(dataMatrix, 2)); % mean-subtract each SVD
% N = size(alpha_trace_mean,2);
% xdft = fft(alpha_trace_mean');
% xdft1 = xdft(1:N/2+1,:);
% nf = size(xdft1,1);
% psdx = (1/(Fs*N)) * abs(xdft1).^2;
% psdx(2:end-1,:) = 2*psdx(2:end-1,:);
% freq1 = 0:Fs/N:Fs/2;
% psdx_mean = mean(psdx,2);
% psdx_std = std(psdx,[],2);
% freq_value = [0.2,0.5,1,2,4,6,8,10,16,32];
% log_freq_value = log10(freq_value);
% freqN = length(psdx_mean);
% 
% plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
% xlim([log_freq_value(1),log_freq_value(end)]);
% xticks(log_freq_value(1:end));
% xticklabels({'0.2','0.5','1','2','4','6','8','10','16'});
% xlabel('log10(Frequency)');
% ylabel('log10(Power) (\muV)^2)');

% hold on;
% dataMatrix = zeros(batchN, epochSize);
% for i = 1:batchN
%     dataMatrix2(i,:) =  meanTrace_filt1(1+epochSize*(i-1):epochSize*i);
% end   
% alpha_trace_mean2 = bsxfun(@minus, dataMatrix2, mean(dataMatrix2, 2)); % mean-subtract each SVD
% N = size(alpha_trace_mean2,2);
% xdft2 = fft(alpha_trace_mean2');
% xdft3 = xdft2(1:N/2+1,:);
% nf = size(xdft3,1);
% psdx2 = (1/(Fs*N)) * abs(xdft3).^2;
% psdx2(2:end-1,:) = 2*psdx2(2:end-1,:);
% freq1 = 0:Fs/N:Fs/2;
% psdx_mean2 = mean(psdx2,2);
% psdx_std2 = std(psdx2,[],2);
% freq_value = [0.2,0.5,1,2,4,6,8,10,16,32];
% log_freq_value = log10(freq_value);
% freqN = length(psdx_mean2);
% 
% plot(log10(freq1(2:freqN)),log10(psdx_mean2(2:freqN)),'color','r');
% xlim([log_freq_value(1),log_freq_value(end)]);
% xticks(log_freq_value(1:end));
% xticklabels({'0.2','0.5','1','2','4','6','8','10','16'});
% xlabel('log10(Frequency)');
% ylabel('log10(Power) (\muV)^2)');

hold on;
dataMatrix = zeros(batchN, epochSize);
for i = 1:batchN
    dataMatrix3(i,:) =  meanTrace_V(1+epochSize*(i-1):epochSize*i);
end   
alpha_trace_mean3 = bsxfun(@minus, dataMatrix3, mean(dataMatrix3, 2)); % mean-subtract each SVD
N = size(alpha_trace_mean3,2);
xdft4 = fft(alpha_trace_mean3');
xdft5 = xdft4(1:N/2+1,:);
nf = size(xdft5,1);
psdx3 = (1/(Fs*N)) * abs(xdft5).^2;
psdx3(2:end-1,:) = 2*psdx3(2:end-1,:);
freq1 = 0:Fs/N:Fs/2;
psdx_mean3 = mean(psdx3,2);
psdx_std3 = std(psdx3,[],2);
freq_value = [0.2,0.5,1,2,4,6,8,10,16,32];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean3);

plot(log10(freq1(2:freqN)),log10(psdx_mean3(2:freqN)),'color','b');
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8','10','16'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (\muV)^2)');
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
imagesc(squeeze(tracePhase1(:,:,frame)));
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
ft = exp(i*(tracePhaseB-tracePhaseA));
%%
ft2 = angle(mean(ft));
%%
figure;
scatter(1:length(tracePhaseA),tracePhaseA,'k');
hold on;
scatter(1:length(tracePhaseA),tracePhaseB,'r');
hold on;
scatter(1:length(tracePhaseA),ft,'b');
%%
frt = mean(ft*35/(2*pi));