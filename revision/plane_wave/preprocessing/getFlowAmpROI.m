function [angle_mean,traceAmp_mean] = getFlowAmpROI(traceAmpt,vxRawt,vyRawt,bwroi)
%%
for frame = 1:size(traceAmpt,3)
    traceAmpta = squeeze(traceAmpt(:,:,frame));
    traceAmptb = traceAmpta(bwroi);
    traceAmp_mean(frame,1) = mean(traceAmptb(:));
end
%% mean flow angle in rectangle roi
vxRaw3_mean = zeros(size(vxRawt,3),1);
vyRaw3_mean = zeros(size(vyRawt,3),1);
for frame = 1:size(vxRawt,3)
    vxRaw1a = squeeze(vxRawt(:,:,frame));
    vyRaw1a = squeeze(vyRawt(:,:,frame));
    vxRaw3 = vxRaw1a(bwroi);
    vyRaw3 = vyRaw1a(bwroi);
    vxRaw3_mean(frame,1) = mean(vxRaw3);
    vyRaw3_mean(frame,1) = mean(vyRaw3);
end
vxy_mean = complex(vxRaw3_mean,vyRaw3_mean);
angle_mean = angle(vxy_mean);