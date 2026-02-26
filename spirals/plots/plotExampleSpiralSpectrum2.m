function [hs1g2] = plotExampleSpiralSpectrum2(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% left SSp 
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
hemi = 'left';
[indexSSp,UselectedSSp] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
%%
kk = 11;                                                                   % LK_0003
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform;
load(fullfile(data_folder,'spirals','spirals_example',[fname '_mask.mat']));
%%
scale = 1;
hist_bin = 40; % pixels
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed = Utransformed(1:8:end,1:8:end,:);
mimgtransformed = mimgtransformed(1:8:end,1:8:end);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate2 = 1;
first_frame =34;
frameStart = 70885-first_frame; frameEnd = frameStart +70;
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
Ur = reshape(Utransformed, size(Utransformed,1)*size(Utransformed,2), size(Utransformed,3));
V1 = V(:,frameTemp);
rawTrace2 = Ur*V1;
rawTrace2 = rawTrace2 -mean(rawTrace2 ,2);
tsize = size(V1,2);
rawTrace2 = reshape(rawTrace2,size(Utransformed,1),size(Utransformed,2),[]);
rawTrace2 = rawTrace2(:,:,1+35/rate2:end-35/rate2);
rawTrace2 = rawTrace2./mimgtransformed;
trace_raw3 = reshape(rawTrace2,size(rawTrace2,1)*size(rawTrace2,2), size(rawTrace2,3));
%%
trace_raw4 = trace_raw3(logical(BW_SSp(:)),:);
index1 = not(isnan(trace_raw4(:,1)));
trace_raw4 = trace_raw4(index1,:);
%%
% clear trace_raw2
% for i = 1:7
%     trace_raw2(i,:) = rawTrace2(pixel(i,1),pixel(i,2),:);
% end
%%
[freq1,psdx,~] = fft_spectrum(trace_raw4);
psdx_mean = mean(psdx,2);
psdx_std = std(psdx,[],2);
freq_value = [0.2,0.5,1,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean);
hs1g2 = figure('Renderer', 'painters', 'Position', [100 100 300 600]);
subplot(1,1,1);
plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
% shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),...
%     log10(psdx_std(2:freqN))./sqrt(size(psdx,2)),'lineProps','k');
hold on;
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8','10'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (df/f^2)');
ylim([-7,-3]);   
print(hs1g2,fullfile(save_folder, 'FigS1g2_example_time_series_spectrum.pdf'),...
    '-dpdf', '-bestfit', '-painters');