function [h1,meanTrace2,tracePhase] = plotExampleTracePhase6(wf_mean2,BW2,maskPath,st,atlas1)
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/2;
%%
h1 = figure('Renderer', 'painters', 'Position', [50 50 950 400]);
for kk = 1:3
    wf_mean = squeeze(wf_mean2(:,:,:,kk));
    wf1 = reshape(wf_mean,size(wf_mean,1)*size(wf_mean,2),size(wf_mean,3));
    meanTrace = wf1 -mean(wf1 ,2);
    meanTrace = double(meanTrace)';
    [f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
    meanTrace = filtfilt(f1,f2,meanTrace);
    meanTrace2 = meanTrace';
    meanTrace2 = reshape(meanTrace2,size(wf_mean,1),size(wf_mean,2),size(wf_mean,3));

    traceHilbert =hilbert(meanTrace);
    tracePhase = angle(traceHilbert);
    tracePhase = tracePhase';
    tracePhase = reshape(tracePhase, size(wf_mean,1),size(wf_mean,2),size(wf_mean,3));
    % mina = min(wf_mean(:));
    % maxa = max(wf_mean(:));
    cmax = 0.03; cmin = -0.03;
    for i = 1:19
        ax1 = subplottight(10,19,(kk-1)*19+i);
        im_raw = imagesc(squeeze(wf_mean(:,:,66+i)));
        % caxis([0.8*mina,0.8*maxa]);
        caxis([cmin,cmax]);
        colormap(ax1,parula);
        axis off; 
        axis image;
        set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);   
        if i==19
            cb = colorbar;
            a =  cb.Position; %gets the positon and size of the color bar
            set(cb,'Position',[a(1) a(2) a(3) a(4)])% To change size
        end

        ax2 = subplottight(10,19,(kk+3)*19+i);
        im_phase = imagesc(tracePhase(:,:,66+i));
        colormap(ax2,colorcet('C06'));
        axis image; axis off;
        set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        caxis([-pi,pi]);
    end
end
% cb4 = subplottight(4,19,19);
% % im1 = imagesc(squeeze(trace2d(:,:,19)),'visible','off');
% im1 = imagesc(squeeze(wf_mean(:,:,66+i)));
% hold on;
% plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
% axis image; axis off;
% set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
% colormap(cb4,'parula');
% caxis([cmin,cmax]);
% axis off;
% cb = colorbar;
% a =  cb.Position; %gets the positon and size of the color bar
% set(cb,'Position',[a(1) a(2) a(3) a(4)])% To change size