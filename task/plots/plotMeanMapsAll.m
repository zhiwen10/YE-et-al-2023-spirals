function hs15bdf = plotMeanMapsAll(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
t2 = -2:1/35:2;
downscale = 16;
scale3 = 5/16;
lineColor = 'k'; lineColor1 = 'w';
labels ={'correct','incorrect', 'miss'};
%% plot correct mean maps
data_folder1 = fullfile(data_folder,'task','task_mean_maps');
load(fullfile(data_folder1,'task_mean_maps_all_mice.mat'));
load(fullfile(data_folder,'tables','task_mask_all_mice.mat'));
downscale = 16;
BW2 = BW2(1:downscale:end,1:downscale:end);
cmax = prctile(trace_mean_all(:),99.98); 
cmin = prctile(trace_mean_all(:),0.02);
%%
hs15bdf = figure('Renderer', 'painters', 'Position', [50 50 950 850]);
for kk = 1:3
    trace_mean3 = squeeze(trace_mean_all(:,:,:,kk));
    tracePhase3 = squeeze(tracePhase_all(:,:,:,kk));
    hemi = [];
    for i = 1:18
        iframe = 70+i;
        ax3 = subplottight(6,19,i+2*19*(kk-1));
        frame_temp = squeeze(trace_mean3(:,:,iframe));
        im1 = imagesc(frame_temp);
        caxis([cmin,cmax]);
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        axis image; axis off;
        set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');    
        colormap(ax3,'parula');
        if kk == 1
            title([num2str(round(t2(iframe)*1000)) 'ms']);
        end

        ax4 = subplottight(6,19,i+19+2*19*(kk-1));
        frame_temp = squeeze(tracePhase3(:,:,iframe));
        im1 = imagesc(frame_temp);
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        axis image; axis off;
        caxis([-pi,pi]);
        set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');   
        colormap(ax4,colorcet('C06'));        
    end
    
    if kk ==1
        cb4 = subplottight(6,19,19);
        im1 = imagesc(squeeze(trace_mean3(:,:,19)),'visible','off');
        colormap(cb4,'parula');
        caxis([cmin,cmax]);
        axis off;
        cb = colorbar;
        a =  cb.Position; %gets the positon and size of the color bar
        set(cb,'Position',[a(1)-0.02 a(2) a(3) a(4)/8])% To change size
    end
    
    dim = [.5 .9-(kk-1)*0.34 .1 .1];
    str = labels{kk};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
end
%%
print(hs15bdf, fullfile(save_folder,['FigS15bdf_mean_maps_all']),...
    '-dpdf', '-bestfit', '-painters');