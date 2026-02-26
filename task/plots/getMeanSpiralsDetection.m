function h5b = getMeanSpiralsDetection(data_folder,save_folder)
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
downscale = 2;
BW2 = BW2(1:downscale:end,1:downscale:end);
%% load mean maps
load(fullfile(data_folder,'task','task_mean_maps',...
    'task_mean_maps_all_mice.mat'));
%%
hemi = [];
scale3 = 5/2;
lineColor = 'k'; lineColor1 = 'w';
v = [36 33; 62 33;62 59;36 59]*8;
f = [1,2,3,4];
%%
kk = 1;
% -2:2 seconds 141 samples
wf_mean2 = squeeze(trace_mean_all(:,:,:,kk));
wf_mean3 = imresize(wf_mean2,[660,570]);
%%
[traceFilt2,tracePhase2] = trace_filt_nan(wf_mean3);
%% spiral detection
freq = [2,8];                                                              % data filtering frequency range
rate = 1;                                                                  % set to 1, if no upsampling in time
params.downscale = 1;                                                      % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0;                                                        % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35;                                                            % frame sampling rate
params.halfpadding = 120;                                                  % half padding should be the max spiral detection radius 
params.padding = 2*params.halfpadding;
params.th = 1:36:360;                                                      % 10 evenly-spaced angles in a circle
params.rs = 10:5:20;                                                       % candidate spiral check at small radii [10, 15, 20 pixels, 17.3um/pixel]
params.gridsize = 10;                                                      % search grid resolution at 10 pixels/grid
params.spiralRange = linspace(-pi,pi,5);                                   % evenly-spaced 4 phase quadrants of a circle
params.gsmooth = 0;                                                        % no spatial smoothing 
params.epochL = 141;                                                       % frames per processing batch 
params.nt = size(tracePhase2,3);                                            % total frames
params.frameRange = 1:params.nt;                                           % frameRange excludeing 1 second edges
params.dThreshold = 15;                                                    % distance threshold for sprial center clusting
params.rsRCheck = 10:10:100;                                               % final spiral radius check range
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);        % roi grids for refined search
%% generate coarse search grids with zeros padded at the edges
params.xsize = size(tracePhase2,1);
params.ysize = size(tracePhase2,2);
xsizePadded = params.xsize+params.padding; 
ysizePadded = params.ysize+params.padding;
[xx,yy] = meshgrid(...
    min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
    min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
params.xx = xx; params.yy = yy;
%% draw brain mask roi
% mimg1 = squeeze(tracePhase2(:,:,100));
% mimg2 = padZeros(mimg1,params.halfpadding);
% figure; 
% ax1 = imagesc(mimg2);
% colormap(colorcet('C06'));
% % roi = drawpolygon;
% % fname1 = 'task_evoked_map_roi.mat';
% % save(fname1,'roi');  
% % draw brain mask roi, or import from roi folder and detect sprials                                           
% %  main search algorithm 
% load(fname1);
% tf = inROI(roi,params.xx(:),params.yy(:));
% params.xxRoi = params.xx(tf);                                                % only use the grids that inside the roi to save time
% params.yyRoi = params.yy(tf);                                                % only use the grids that inside the roi to save time         
% %% spiral detection
% pwAll = [];
% pwAll1 = []; pwAll2 = []; 
% pwAll3 = []; pwAll4 = []; pwAll5 = [];
% tracePhase1 = padZeros(tracePhase2,params.halfpadding);                       % pad tracephase with edge zeros
% nframe = size(tracePhase1,3);
% frameStart = 1;
% for frame = 1:nframe
%     A = squeeze(tracePhase1(:,:,frame));
%     pwAll1 = spiralAlgorithm(A,params);                                      % coarse search
%     frameN = frame+frameStart-1;
%     pwAll2 = checkClusterXY(pwAll1,params.dThreshold);                      % spatial clustering of nearby duplicate spiral centers
%     pwAll3 = doubleCheckSpiralsAlgorithm(A,pwAll2,params);                  % double check the mean cluster centers from the last step are still spiral centers   
%     pwAll4 = spatialRefine(A,pwAll3,params);                                % refined spiral search based on candidate sprial centers
%     [pwAll5] = spiralRadiusCheck2(A,pwAll4,params);                         % find spiral radius and traveling direction
%     if not(isempty(pwAll5))
%         pwAll5(:,end+1) = frameN;                                            % attach frame ID label
%     end
%     pwAll = [pwAll;pwAll5];                                                  % concatenate all sprials from the current batch
% end
% pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;                            % recalculate spiral 2d coordinates without padding   
% 
% fname = 'task_spirals.csv';                                             % save data file name
% T1 = array2table(pwAll,'VariableNames',...                                 % make the sprial matrix into a table with varibale names
%     {'spiral_center_x','spiral_center_y',...
%     'spiral_radius','spiral_direction','spiral_frame'});            
% writetable(T1,fname);                            % write to csv file
% %% temporal grouping
% T2 = readtable(fname);
% pwAll = table2array(T2);
% filteredSpirals = pwAll(pwAll(:,3)>=40,:);                                   % only use sprials with radius >40 pixels, based on 3d-fft
% filteredSpirals =unique(filteredSpirals, 'rows');                            % get rid of duplication, in case any
% filteredSpirals = sortrows(filteredSpirals,5);                               % sort based on frame number
% [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);            % main grouping algorithm
% save('task_spirals_group_fftn.mat','archiveCell');
%%
% cmax = prctile(wf_mean3(:),99.98); 
% cmin = prctile(wf_mean3(:),0.02);
cmin = -0.01;
cmax = 0.01;
%%
load('task_spirals_group_fftn.mat','archiveCell');
spirals = cell2mat(archiveCell);
h5b = figure('Renderer', 'painters', 'Position', [50 50 950 400]);
th2 = 1:5:360; 
for i = 1:14
    iframe = i+68;
    spiral_temp = spirals(spirals(:,5) == 68+i,:);
    ax1 = subplottight(4,15,i);
    frame_temp1 = squeeze(wf_mean3(:,:,iframe));
    im1 = imagesc(frame_temp1);
    caxis([cmin,cmax]);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');    
    colormap(ax1,'parula');
    hold on;
    patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);

    ax2 = subplottight(4,15,i+15);
    frame_temp2 = squeeze(tracePhase2(:,:,iframe));
    im1 = imagesc(frame_temp2);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    caxis([-pi,pi]);
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');   
    colormap(ax2,colorcet('C06'));     
    hold on;
    patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);

    %%
    ax3(i) = subplottight(4,15,i+30);
    ax3(i).Position(1) = ax3(i).Position(1);
    ax3(i).Position(2) = ax3(i).Position(2);
    ax3(i).Position(3) = ax3(i).Position(3)-0.01;
    ax3(i).Position(4) = ax3(i).Position(4)-0.01;
    im1 = imagesc(frame_temp1);
    colormap(ax3(i),'parula'); 
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    caxis([cmin,cmax]);
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    xlim([36,62]*8);
    ylim([33,59]*8);    
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end 
    %%
    ax4(i) = subplottight(4,15,i+45);
    ax4(i).Position(1) = ax4(i).Position(1);
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.01;
    ax4(i).Position(4) = ax4(i).Position(4)-0.01;
    im1 = imagesc(frame_temp2);
    colormap(ax4(i),colorcet('C06'));  
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    caxis([-pi,pi]);
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    xlim([36,62]*8);
    ylim([33,59]*8);     
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end 
end
cb4 = subplottight(4,15,15);
im1 = imagesc(squeeze(wf_mean3(:,:,82)),'visible','on');
colormap(cb4,'parula');
caxis([cmin,cmax]);
axis off;
cb = colorbar;
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)-0.02 a(2) a(3) a(4)/4])% To change size

%%
print(h5b,fullfile(save_folder,'Fig5b_correct_mean_maps4.pdf'),...
    '-dpdf', '-bestfit', '-painters');