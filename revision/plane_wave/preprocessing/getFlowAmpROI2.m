function [sync,traceAmp_mean] = getFlowAmpROI2(traceAmpt,vxRawt,vyRawt,bwroi)
%%
for frame = 1:size(traceAmpt,3)
    traceAmpta = squeeze(traceAmpt(:,:,frame));
    traceAmptb = traceAmpta(bwroi);
    traceAmp_mean(frame,1) = mean(traceAmptb(:));
end
%%
vxyt = complex(vxRawt,vyRawt);
vxyt2  = reshape(vxyt,size(vxyt,1)*size(vxyt,2),size(vxyt,3));
vxyt3 = vxyt2(logical(bwroi(:)),:);
%%    
vxyt3 = vxyt3./abs(vxyt3);                                                 % normalize to unit vector
sync = sum(vxyt3,1)./size(vxyt3,1);                                        % get sync index