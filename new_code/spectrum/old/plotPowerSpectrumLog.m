function plotPowerSpectrumLog(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
BW = logical(projectedAtlas1);
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
pixels(1,:) = [845,835]; % VISp
pixels(2,:) = [775,650]; % RSP
pixels(3,:) = [590,750]; % SSp-ul
pixels(4,:) = [520,850]; % SSp-ll
pixels(5,:) = [480,950]; % SSp-m
pixels(6,:) = [550,950]; % SSp-n
pixels(7,:) = [675,905]; % SSp-bfd
area_names = {'VISp','RSP','SSp-ul','SSp-ll','SSp-m','SSp-n','SSp-bfd'};
%%
color2 = cbrewer2('seq','YlOrRd',9);
scale3 = 5; hemi = []; lineColor = 'w';
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 700]);
% params
color3 = cbrewer2('qual','Set3',15);
for kk = 1:15
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum2',...
        [fname '_fftSpectrum2.mat']));
    % session_root = fullfile(data_folder,'spirals\svd',fname);
    % mimg = readNPY(fullfile(session_root, 'meanImage.npy'));
    % load(fullfile(data_folder,'spirals\rf_tform',...
    %     [fname '_tform.mat']));                                                % load atlas transformation matrix tform;
    % mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1))); 
    pixels_copy = pixels;
    if kk == 3
        pixels_copy(5,:) = [560,920];
        pixels_copy(6,:) = [655,890];
    end
    xticks_value = log10([0.5,2,4,6,8]);
    for i = 1:7
        ax(i) = subplot(2,4,1+i);
        powerV = log10(squeeze(psdxMeanTransformed(pixels_copy(i,1),pixels_copy(i,2),:)))';
        plot(log10(freq(2:end)),powerV(2:end),'color',color3(kk,:));
        xlabel('Frequency (Hz)'); ylabel('log10(Power)');
        xticks(xticks_value);
        xticklabels({'0.5','2','4','6','8'});
        hold on;  
    end
end
for i = 1:7
    ax(i) = subplot(2,4,1+i);
    title(area_names{i});
end
% get mimg for an example session
kk = 1;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals\svd',fname);
mimg = readNPY(fullfile(session_root, 'meanImage.npy'));
load(fullfile(data_folder,'spirals\rf_tform',...
    [fname '_tform.mat']));                                                % load atlas transformation matrix tform;
mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1))); 

ax0 = subplot(2,4,1);
im_raw = imagesc(mimgt);
set(im_raw, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
axis image; axis off;
hold on;
overlayOutlines(coords,1);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(ax0,pixels(:,2),pixels(:,1),16,color2(3:end,:),'filled');
%%
print(h1, fullfile(save_folder,'power_spectrum_all_log'), '-dpdf', '-bestfit', '-painters');
