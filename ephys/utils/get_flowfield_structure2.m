function [data_raw] = get_flowfield_structure2(Ut1,dV1,epochs,mimg,flow)
%%
if nargin <5
    flow = 1;
end
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
t = [0,1/35,2/35];
useGPU = 0;
frameStart = epochs(1); frameEnd = epochs(end);
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV_raw_epoch = dV1(1:50,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(Ut1(:,:,1:50),dV_raw_epoch,t,params,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimg; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
traceAmp1 = traceAmp1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data   

trace2d1 = permute(trace2d1,[3,1,2]);
tracePhase1 = permute(tracePhase1,[3,1,2]);
traceAmp1 = permute(traceAmp1,[3,1,2]);


if flow
    [vxRaw,vyRaw] = HS_flowfield(tracePhase1,useGPU);
    data_raw.vxRaw = vxRaw;
    data_raw.vyRaw = vyRaw;
end

tracePhase1 = tracePhase1(1:end-1,:,:);
traceAmp1 = traceAmp1(1:end-1,:,:);
trace2d1 = trace2d1(1:end-1,:,:);

data_raw.tracePhase1 = tracePhase1;
data_raw.traceAmp1 = traceAmp1;
data_raw.trace2d1 = trace2d1;

end