function h5b = plotCorrectMapsFlow2(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%%
load(fullfile(data_folder,'tables','task_mask_all_mice.mat'));
downscale = 16;
BW2 = BW2(1:downscale:end,1:downscale:end);
%% load mean maps
load(fullfile(data_folder,'task','task_mean_maps',...
    'task_mean_maps_all_mice.mat'));
%%
hemi = [];
scale3 = 5/16;
lineColor = 'k'; lineColor1 = 'w';
v = [37 38; 58 38;58 59;37 59];
f = [1,2,3,4];
%
h5b = figure('Renderer', 'painters', 'Position', [50 50 950 600]);
for kk = 1
    % -2:2 seconds 141 samples
    trace_mean_current = squeeze(trace_mean_all(:,:,:,kk));
    %%
    [traceFilt2,tracePhase2] = trace_filt_nan(trace_mean_current);
    tracePhase3 = permute(tracePhase2,[3,1,2]);
    tracePhase3 = tracePhase3(:,35:59,10:60);
    useGPU = 0;
    [vxRaw_all,vyRaw_all] = HS_flowfield(tracePhase3,useGPU);
    %% only use -1:1 seconds 70 samples
    trace_mean4 = trace_mean_current(:,:,54:88); 
    traceFilt4 = traceFilt2(:,:,54:88); 
    tracePhase4 = tracePhase2(:,:,54:88);
    vxRaw_all1 = vxRaw_all(54:88,:,:);
    vyRaw_all1 = vyRaw_all(54:88,:,:);
    %%
    if kk == 1
        cmax = prctile(trace_mean4(:),99.98); 
        cmin = prctile(trace_mean4(:),0.02);
    end
    %%
    for i = 1:14
        iframe = i+15;
        ax1 = subplottight(3,15,(kk-1)*15*3+i);
        frame_temp1 = squeeze(trace_mean4(:,:,iframe));
        im1 = imagesc(frame_temp1);
        caxis([cmin,cmax]);
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        axis image; axis off;
        set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');    
        colormap(ax1,'parula');

        ax2 = subplottight(3,15,(kk-1)*15*3+i+15);
        frame_temp2 = squeeze(tracePhase4(:,:,iframe));
        im1 = imagesc(frame_temp2);
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        axis image; axis off;
        caxis([-pi,pi]);
        set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');   
        colormap(ax2,colorcet('C06'));     
        hold on;
        patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);
    end
    cb4 = subplottight(3,15,15);
    im1 = imagesc(squeeze(trace_mean4(:,:,18)),'visible','off');
    colormap(cb4,'parula');
    caxis([cmin,cmax]);
    axis off;
    cb = colorbar;
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)-0.02 a(2) a(3) a(4)/8])% To change size
end
%%
print(h5b,fullfile(save_folder,'Fig5b_correct_mean_maps2.pdf'),...
    '-dpdf', '-bestfit', '-painters');