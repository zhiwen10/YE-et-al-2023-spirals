function [tracePhase1_raw1,vxy_raw1,traceAmp_raw1] =get_flowfield5(Ut,dV_raw,mimg1,flow,len)
%% angle difference for oscillation period  
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
t1 = [0,1/35,2/35];
useGPU = 1;
%%
frameN = size(dV_raw,2);
if frameN > len
    frameN = len;
end
%%
dV_raw1 = squeeze(dV_raw(:,1:frameN));
[~,traceAmp,tracePhase1_raw1] = spiralPhaseMap4(Ut,dV_raw1,t1,params,rate,mimg1);
traceAmp_raw1 = permute(traceAmp,[3,1,2]);
tracePhase1_raw1 = permute(tracePhase1_raw1,[3,1,2]);
%%
if flow
    [vxRaw,vyRaw] = HS_flowfield(tracePhase1_raw1,useGPU);
else
    vxRaw = []; vyRaw = [];
    vxPred = []; vyPred = [];
end
%%
vxy_raw1 = complex(vxRaw, vyRaw);
traceAmp_raw1 = traceAmp_raw1(1:end-1,:,:);
tracePhase1_raw1 = tracePhase1_raw1(1:end-1,:,:);
