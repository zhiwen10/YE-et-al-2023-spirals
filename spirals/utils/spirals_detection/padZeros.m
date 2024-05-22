function [tracePhase] = padZeros(tracePhase1,halfpadding)
% tracePhase size should be xsize * ysize * nframes
padding = 2* halfpadding;
xsize = size(tracePhase1,1); 
ysize = size(tracePhase1,2);
nframe = size(tracePhase1,3);
tracePhase = zeros(xsize+padding ,ysize+padding,nframe);
tracePhase(1+halfpadding:halfpadding+xsize,1+halfpadding:halfpadding+ysize,:) = tracePhase1;
end