function hs13a = plotVarMap(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%% cortex surface outline
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
prediction_folder = fullfile(data_folder,'ephys','dv_prediction');
roi_folder = fullfile(data_folder,'ephys','roi');
regist_folder = fullfile(data_folder,'ephys','rf_tform_4x');
%%
scale = 1;
scale3 = 5/scale;
hemi = [];
hs13a = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);
count1 = 1;
for areai = [1,2,4]
    current_T = T(contains(T.Area,area(1,areai)),:);
    for kk = 1:size(current_T)
        pos = kk+(count1-1)*12;
        ax2 = subplottight(3,12,pos);
        ops = get_session_info2(current_T,kk,data_folder);
        fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
        fname_predict = fullfile(prediction_folder,[fname '_dv_predict.mat']);
        load(fname_predict,'explained_var_all');
        load(fullfile(roi_folder,[fname '_roi.mat']));
        explained_var_all(explained_var_all<0) = 0;
        var_mean = squeeze(mean(explained_var_all,3));
        fname1 = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x'];
        load(fullfile(regist_folder, fname1));
        mean_var_t = imwarp(var_mean,tform,'OutputView',imref2d(size(projectedTemplate1)));
        im2 = imagesc(mean_var_t);
        set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
        caxis([0,0.8]);
        axis image; axis off;
        cmap1 = inferno;
        colormap(ax2,cmap1);
        
        hold on;
        plotOutline(maskPath(1:3),st,atlas1,hemi,scale3);
        plotOutline(maskPath(4),st,atlas1,[],scale3);
        plotOutline(maskPath(5),st,atlas1,[],scale3);
        plotOutline(maskPath(6:11),st,atlas1,[],scale3);
        text(0,50,ops.mn,'Interpreter', 'none');
    end
    count1 = count1+1;
end
ax3 = subplottight(3,12,35);
hold off;
im2 = imagesc(mean_var_t);
set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
caxis([0,0.8]);
axis image; axis off;
cmap1 = inferno;
colormap(ax3,cmap1);
hold on;
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
cb = colorbar;
cb.Ticks = [0,0.4,0.8];
cb.TickLabels = string([0,0.4,0.8]);
%%
print(hs13a, fullfile(save_folder,'FigS13a_all_variance_maps.pdf'),...
    '-dpdf', '-bestfit', '-painters');