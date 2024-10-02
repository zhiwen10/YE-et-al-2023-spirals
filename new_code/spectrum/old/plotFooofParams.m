function plotFooofParams(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
BW2 = imresize(BW,1/8);
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%% plot fooof params
data_folder1 = fullfile(data_folder,'spirals\spectrum\fooof_sessions');
session{1} = 1:8;
session{2} = 9:15;
for sessioni = 1:2
    h1 = figure('Renderer', 'painters', 'Position', [10 10 900 1000]);
    count1 = 1;
    for kk = session{sessioni}                                                             
        mn = T.MouseID{kk};
        tda = T.date(kk);
        en = T.folder(kk);    
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        fname = [mn '_' tdb '_' num2str(en)];
        %%
        Ta = readtable(fullfile(data_folder1,[fname '_fooof.csv']));
        fooofV1 = table2array(Ta);
        indx = readtable(fullfile(data_folder1,[fname '_index.csv']));
        indx1 = logical(indx.Index);
        %%
        fooofV = nan(165*143,8);
        fooofV(indx1,:) = fooofV1;
        fooofV = reshape(fooofV,143,165,[]);
        %%
        scale3 = 5/8; lineColor = 'w'; hemi = [];
        ax1 = subplot(8,3,1+(count1-1)*3);
        exponents = squeeze(1./fooofV(:,:,3))';
        exponents(exponents>1) = 1;
        exponents(exponents<0) = 0;
        im1 = imagesc(exponents);
        colormap(ax1, parula);
        axis image; axis off;
        set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        colorbar;
        title('Exponents');
        text(-50,-50,mn,'Interpreter','None');
        
        ax2 = subplot(8,3,2+(count1-1)*3);
        im2 = imagesc(squeeze(fooofV(:,:,4))');
        colormap(ax2, parula);
        axis image; axis off;
        set(im2, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        colorbar;
        title('Central Frequency (Hz)');
        ax3 = subplot(8,3,3+(count1-1)*3);
        im3 = imagesc(squeeze(fooofV(:,:,5))');
        colormap(ax3, parula);
        axis image; axis off;
        set(im3, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        colorbar;
        title('Peak Power');
        %%
        count1 = count1+1;
    end
    print(h1, fullfile(save_folder,['fooof_' num2str(sessioni) '.pdf']),...
    '-dpdf', '-bestfit', '-painters');
end