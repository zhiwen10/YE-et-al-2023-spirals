%%
folder = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git\YE-et-al-2023-spirals_20260220\simulation\files2\phase_data';
phasemap1 = readNPY(fullfile(folder,'phase_evolution_cartesian.npy'));
phasemap2 = readNPY(fullfile(folder,'phase_evolution_polar.npy'));
%%
seed = 3;
phasemapA = squeeze(phasemap1(seed,:,:,:));
phasemapB = squeeze(phasemap2(seed,:,:,:));
h1b = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);

for i = 1:10
    ax1 = subplottight(2,10,i);
    ax1.Position(1) = ax1.Position(1)+0.005;
    ax1.Position(2) = ax1.Position(2);
    ax1.Position(3) = ax1.Position(3)-0.01;
    ax1.Position(4) = ax1.Position(4)-0.01;
    imagesc(squeeze(phasemapA(i*20,:,:)));
    colormap(ax1,colorcet('C06'));
    axis image; axis off;
    
    ax2 = subplottight(2,10,10+i);
    ax2.Position(1) = ax2.Position(1)+0.005;
    ax2.Position(2) = ax2.Position(2);
    ax2.Position(3) = ax2.Position(3)-0.01;
    ax2.Position(4) = ax2.Position(4)-0.01;
    imagesc(squeeze(phasemapB(i*20,:,:)));
    colormap(ax2,colorcet('C06'));
    axis image; axis off;
    
end
%%
print(h1b, 'spiral_circular_example.pdf','-dpdf', '-bestfit', '-painters');
%%
pwAlla = cell(7,1);
for seed = 1:7
    %%
    % phasemapA = squeeze(phasemap1(seed,:,:,:));
    phasemapA = squeeze(phasemap2(seed,:,:,:));
    phasemapA = permute(phasemapA,[2,3,1]);
    
    params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
    params.lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
    params.Fs = 35; % frame sampling rate
    params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
    params.padding = 2*params.halfpadding;
    params.th = 1:36:360; 
    params.rs = 10:5:20;
    params.gridsize = 10;
    params.spiralRange = linspace(-pi,pi,5);
    params.gsmooth = 0;
    frameN1 = size(phasemapA,3);
    params.dThreshold = 15;
    params.rsRCheck = 10:10:100;
    params.xsize = size(phasemapA,1);
    params.ysize = size(phasemapA,2);
    [params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10); 
    xsizePadded = params.xsize+params.padding; 
    ysizePadded = params.ysize+params.padding;
    [xx,yy] = meshgrid(min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1,...
        min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
    params.xx = xx; params.yy = yy;

    phasemapA1 = phasemapA;
    phasemapA1 = padZeros(phasemapA1,params.halfpadding); 
    %% draw brain mask roi
%     mimg1 = squeeze(phasemapA(1:params.downscale:end,1:params.downscale:end,1));
%     mimg2 = padZeros(mimg1,params.halfpadding);
%     figure; 
%     ax1 = imagesc(mimg2);
%     roi = drawpolygon;
%     fname1 = 'model_roi.mat';
%     save(fname1,'roi');                                            % results saved in <full_roi> folder 
    load(fullfile(folder,'model_roi.mat'));
    tf = inROI(roi,params.xx(:),params.yy(:));
    params.xxRoi = params.xx(tf);                                          % only use the grids that inside the roi to save time
    params.yyRoi = params.yy(tf); 
    %%
    pwAll = []; pwAll1 = []; pwAll2 = []; pwAll3 = []; pwAll4 = [];
    tic
    for frame = 1:frameN1
        A = squeeze(phasemapA1(:,:,frame));
        A1 = wrapToPi(A(:));
        A = reshape(A1,size(A,1),size(A,2));
        [pwAll1] = spiralAlgorithm(A,params);
        pwAll2 = checkClusterXY(pwAll1,params.dThreshold);
        [pwAll3] = doubleCheckSpiralsAlgorithm(A,pwAll2,params);
        [pwAll4] = spiralRadiusCheck2(A,pwAll3,params);
        if not(isempty(pwAll4))
            pwAll4(:,end+1) = frame; 
        end
        pwAll = [pwAll;pwAll4];
        if mod(frame,50) == 0
            fprintf('Frame %g/%g; time elapsed %g seconds \n', [frame,frameN1, toc])
        end
    end
    if not(isempty(pwAll))
        pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;
    end
    %%
    pwAlla{seed,1} = pwAll;
end
save('circular_spirals' ,'pwAlla');
%%
folder = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git\YE-et-al-2023-spirals\simulation\spirals';
%%
th2 = 1:5:360;
r = 256;
px1 = 256; py1 = 256;
cx2 = round(r*cosd(th2)+px1);
cy2 = round(r*sind(th2)+py1);

color1 = cbrewer2('qual','Set1',7);
h1b = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
load(fullfile(folder,'isotropic_spirals.mat'));
for k = 1:5
    ax1 = subplot(2,5,k);
    for i = 1:7
        clear pwAll pwAll2
        pwAll = pwAlla{i};
        pwAll2 = pwAll((pwAll(:,5)>=20*(k-1) & pwAll(:,5)<20*k) ,:);
        scatter(pwAll2(:,1),pwAll2(:,2),8,'markerFaceColor','None','markerEdgeColor',color1(i,:));
        hold on;
    end
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k-','LineWidth',2); 
    axis equal;
    xlim([0,512]);
    ylim([0,512]);
    axis off;
end


load(fullfile(folder,'circular_spirals.mat'));
for k = 1:5
    ax2 = subplot(2,5,k+5);
    for i = 1:7
        clear pwAll pwAll2
        pwAll = pwAlla{i};
        pwAll2 = pwAll((pwAll(:,5)>=20*(k-1) & pwAll(:,5)<20*k) ,:);
        scatter(pwAll2(:,1),pwAll2(:,2),8,'markerFaceColor','None','markerEdgeColor',color1(i,:));
        hold on;
    end
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k-','LineWidth',2); 
    axis equal;
    xlim([0,512]);
    ylim([0,512]);
    axis off;
end
%%
print(h1b, 'circular_model_centers.pdf','-dpdf', '-bestfit', '-painters');