function hr5fg = plotExampleAlphaEpoch(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
pixel(1,:) = [845,835]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,950]; % SSp-m
pixel(6,:) = [550,950]; % SSp-n
pixel(7,:) = [682,905]; % SSp-bfd
pixel(8,:) = [290,700]; % MOs
color1 = cbrewer2('qual','Set3',8);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd','MOs'};
%%
kk = 7; % ZYE_0012
% session info
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(td,'yyyymmdd');
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
dV = dV(1:50,:);
% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,1:50);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
%
pixel_copy = pixel;
if kk == 3
    pixel_copy(5,:) = [560,920];
    pixel_copy(6,:) = [655,890];
end
pixel_copy = round(pixel_copy/params.downscale);
%%    
for i = 1:8
    rawTrace(i,:) = squeeze(Utransformed(pixel_copy(i,1),pixel_copy(i,2),:))'*dV(1:50,:);
    rawTrace(i,:) = rawTrace(i,:)./mimgtransformed(pixel_copy(i,1),pixel_copy(i,2));
end
%
 for i = 1:8
    rawTraceV(i,:) = squeeze(Utransformed(pixel_copy(i,1),pixel_copy(i,2),:))'*V(1:50,:);
    rawTraceV(i,:) = rawTraceV(i,:)./mimgtransformed(pixel_copy(i,1),pixel_copy(i,2));
end
%%
Fs = 35;
rawTrace = double(rawTrace);
rawTrace = rawTrace -mean(rawTrace ,2); 
freq = [2,8];
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
traceFilt = filtfilt(f1,f2,rawTrace');
traceFilt = traceFilt';
traceHilbert =hilbert(traceFilt);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);    
traceAmp2 = movmean(traceAmp,17,2); % smooth mean of 0.5s
traceAmp3 = traceAmp2(3,:);
traceFilt3 = traceFilt(3,:);
%%    
threshold = 0.0025;
[alpha_epoch2,alpha_binary] = getAlphaEpoch(traceAmp3,threshold);
traceFilt3b = traceFilt3;
traceFilt3b(not(alpha_binary)) = nan;
alpha_ratio = sum(alpha_binary)./numel(alpha_binary);
%
hr5fg = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,3,[1,2]);
plot(t,traceFilt3,'k');
hold on;
plot(t,traceAmp3,'r');
hold on;
yline(threshold,'g--','lineWidth',2);
hold on;
plot(t,traceFilt3b,'b');
xlim([1955,1980]);
xticks([1955:5:1980])
xticklabels({'0','5','10','15','20','25'});
subplot(1,3,3);
histogram(traceAmp,'DisplayStyle','stairs');
xline(threshold);
%%
print(hr5fg, fullfile(save_folder,'FigR5fg_exampleAlphaEpoch'),...
    '-dpdf', '-bestfit', '-painters');