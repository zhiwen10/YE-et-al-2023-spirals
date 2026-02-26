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
kk = 6;
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
meanTrace_v = Ur*V(1:50,1:60000);
meanTrace_v = meanTrace_v -mean(meanTrace_v ,2);
meanTrace_v = double(meanTrace_v)';
traceHilbert_v =hilbert(meanTrace_v);
tracePhase_v = angle(traceHilbert_v);
clear traceHilbertV
tracePhase_v = reshape(tracePhase_v,[],size(U2,1),size(U2,2));
tracePhase_v = permute(tracePhase_v,[2,3,1]);
meanTrace_v = reshape(meanTrace_v,[],size(U2,1),size(U2,2));
%% raw phase dv
meanTrace_dv = Ur*dV(1:50,1:60000);
meanTrace_dv = meanTrace_dv -mean(meanTrace_dv ,2);
meanTrace_dv = double(meanTrace_dv)';
traceHilbert_dv =hilbert(meanTrace_dv);
tracePhase_dv = angle(traceHilbert_dv);
clear traceHilbert_dv
tracePhase_dv = reshape(tracePhase_dv,[],size(U2,1),size(U2,2));
tracePhase_dv = permute(tracePhase_dv,[2,3,1]);
meanTrace_dv = reshape(meanTrace_dv,[],size(U2,1),size(U2,2));
%% filtered phase
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
meanTrace_filt = filtfilt(f1,f2,meanTrace_dv);
traceHilbert_filt =hilbert(meanTrace_filt);
tracePhase_filt = angle(traceHilbert_filt);
clear traceHilbert_filt
tracePhase_filt = reshape(tracePhase_filt,[],size(U2,1),size(U2,2));
tracePhase_filt = permute(tracePhase_filt,[2,3,1]);
meanTrace_filt = reshape(meanTrace_filt,[],size(U2,1),size(U2,2));
%%
figure;
imagesc(mimg(200:4:400,250:4:500));
%%
figure;
plot(t(1:60000),squeeze(meanTrace_dv(:,35,30)))
%%
pixel = [200+35*4,250+30*4];
% epoch = [1005,1020];
% epoch = [195,210];
epoch = [1008,1018];
frameStart = find(t>epoch(1),1,'first'); frameEnd = find(t>epoch(2),1,'first');
sample_epoch = frameStart:frameEnd;
%%
sample_epoch2 = sample_epoch;
% sample_epoch2 = sample_epoch-sample_epoch(1)+1;
%% trace example
h1b = figure('Renderer', 'painters', 'Position', [100 100 500 800]);
subplot(6,1,1);
plot(t(sample_epoch2)-t(sample_epoch2(1)),squeeze(meanTrace_v(sample_epoch,35,30)));
xticks([0:2:10]);
xlim([0,10]);
box off;

subplot(6,1,2);
tracePhase_v1 = squeeze(tracePhase_v(35,30,sample_epoch));
tracePhase_v1 = unwrap(tracePhase_v1);
instFreq_v1 = abs(diff(tracePhase_v1))*35/(2*pi);
instFreqSmooth_v1 = smoothdata(instFreq_v1,"gaussian",17);

plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreq_v1,'r');
hold on;
plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreqSmooth_v1,'k');
ylim([0,3]);
box off;
xticks([0:2:10]);
xlim([0,10]);


subplot(6,1,3);
plot(t(sample_epoch2)-t(sample_epoch2(1)),squeeze(meanTrace_dv(sample_epoch,35,30)));
xticks([0:2:10]);
xlim([0,10]);
box off;

subplot(6,1,4);
tracePhase_dv1 = squeeze(tracePhase_dv(35,30,sample_epoch));
tracePhase_dv1 = unwrap(tracePhase_dv1);
instFreq_dv1 = abs(diff(tracePhase_dv1))*35/(2*pi);
instFreqSmooth_dv1 = smoothdata(instFreq_dv1,"gaussian",17);

plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreq_dv1,'r');
hold on;
plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreqSmooth_dv1,'k');
xlim([0,10]);
xticks([0:2:10]);
box off;

subplot(6,1,5);
plot(t(sample_epoch2)-t(sample_epoch2(1)),squeeze(meanTrace_filt(sample_epoch,35,30)));
xticks([0:2:10]);
xlim([0,10]);
box off;

subplot(6,1,6);
tracePhase_filt1 = squeeze(tracePhase_filt(35,30,sample_epoch));
tracePhase_filt1 = unwrap(tracePhase_filt1);
instFreq_filt1 = abs(diff(tracePhase_filt1))*35/(2*pi);
instFreqSmooth_filt1 = smoothdata(instFreq_filt1,"gaussian",17);

plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreq_filt1,'r');
hold on;
plot(t(sample_epoch2(1:end-1))-t(sample_epoch2(1)),instFreqSmooth_filt1,'k');
xticks([0:2:10]);
xlim([0,10]);
box off;
%%
print(h1b, 'Instantenous frequency_example.pdf','-dpdf', '-bestfit', '-painters');
%%
figure;
histogram(instFreq_filt1)
%%
meanTrace2 = Ur*dV(1:50,:);
meanTrace2 = meanTrace2 -mean(meanTrace2 ,2);
meanTrace2 = double(meanTrace2)';
meanTrace_raw2 = reshape(meanTrace2,[],size(U2,1),size(U2,2));
meanTrace_raw3 = squeeze(meanTrace_raw2(:,35,30));
%%
meanTrace_raw1 = squeeze(meanTrace_dv(:,15,30));
meanTrace_filt1 = squeeze(meanTrace_filt(:,15,30));
Fs = 35;
epochSize = 20*Fs;
batchN = floor(size(meanTrace_raw1,1)/epochSize);
clear dataMatrix3
h1b = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
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
    dataMatrix3(i,:) =  meanTrace_v(1+epochSize*(i-1):epochSize*i);
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
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
frameTemp = [35170:35250];
% frameTemp = [35200:35500];
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U,dV1,t,params,freq,rate);
%%
% frame_i = 35297-35200+1;
frame_i = 35211-35170+1;
% frame_i = 35400-35200+1;
figure;
for k = 1:10
    ax1 = subplottight(1,10,k);
    imagesc(squeeze(tracePhase1(:,:,frame_i+k)));
    colormap(ax1,colorcet('C06'));
    axis image; axis off;
end
%%
th2 = 1:5:360; 
h2 = figure('Renderer', 'painters', 'Position', [100 100 800 400]);
ax1 = subplot(1,2,1);
% imagesc(squeeze(tracePhase1(:,:,frame_i)));
% colormap(ax1,colorcet('C06'));
imagesc(mimg);
colormap(gray);
axis image; axis off;
hold on;
rectangle('position',[250,200,250,200],'EdgeColor','w');
hold on;
scatter(pixel(2),pixel(1),24,'r');
ax2 = subplot(1,2,2);
frame = 35211;
imagesc(squeeze(tracePhase1(:,:,frame_i)));
colormap(ax2,colorcet('C06'));
axis image; axis off;
hold on;
spiral_temp = filteredSpirals(filteredSpirals(:,5) == frame,:);
scale1 = 1;
if not(isempty(spiral_temp))
    for kk = 1:size(spiral_temp,1)
        px1 = spiral_temp(kk,1)/scale1;
        py1 = spiral_temp(kk,2)/scale1;
        r = spiral_temp(kk,3)/scale1;
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
        scatter(spiral_temp(kk,1)/scale1,spiral_temp(kk,2)/scale1,8,color1,'filled');
    end
end 
xlim([250,500]);
ylim([200,400]);
%%
print(h2, 'Instantenous frequency_example_frame.pdf','-dpdf', '-bestfit', '-painters');
%%
cx3 = cx2(cx2<63 &cy2<51);
cy3 = cy2(cx2<63 &cy2<51);
for k = 1:length(cx3)
    tracePhaseA(k,1) = tracePhase_dv(cy3(k),cx3(k),82);
    tracePhaseB(k,1) = tracePhase_dv(cy3(k),cx3(k),83);
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