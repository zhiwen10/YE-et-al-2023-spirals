function params = setSpiralDetectionParamsFFTN(mimg,t,roi_ap,roi_ml)
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
params.epochL = 100;                                                       % frames per processing batch 
params.nt = numel(t);                                                      % total frames
params.frameN1 = round((params.nt-70)/100)*100;                            % round-up frame numbers excluding 1 second at both ends for filtering edges
params.frameRange = 36:params.epochL:(36+params.frameN1-1);                % frameRange excludeing 1 second edges
params.dThreshold = 15;                                                    % distance threshold for sprial center clusting
params.rsRCheck = 10:10:100;                                               % final spiral radius check range
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);        % roi grids for refined search
%% generate coarse search grids with zeros padded at the edges
mimg = mimg(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2));
params.xsize = size(mimg,1);
params.ysize = size(mimg,2);
xsizePadded = params.xsize+params.padding; ysizePadded = params.ysize+params.padding;
[xx,yy] = meshgrid(min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
    min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
params.xx = xx; params.yy = yy;
%% apply mask, this helps speed up spiral detection later
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
mimg2 = padZeros(mimg1,params.halfpadding);
mimg2_bw = logical(mimg2);
tf = zeros(size(xx));
for i = 1:numel(xx(:))
    tf(i) = mimg2_bw(xx(i),yy(i));
end
tf = logical(tf);
params.xxRoi = yy(tf); 
params.yyRoi =xx(tf);