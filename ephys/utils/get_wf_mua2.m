function [Ut,mimg,V1,dV1,MUA_std] =get_wf_mua2(ops)
%% get UVT and sp
[U,V,t,mimg] = loadUVt1(ops.session_root);  
dV = [zeros(size(V,1),1) diff(V,[],2)];
% load spike data
[sp] = loadKSdir2(ops.session_root);
%% load tform and variance map
fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
%% registration for phasemap
scale = 4;
Ut = U(1:scale:end,1:scale:end,1:50);
%% prediction from spiking data
[syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t);
WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
dV1 = double(dV(:,~isnan(WF2ephysT1)));
V1 = double(V(1:50,~isnan(WF2ephysT1)));
[MUA_std] = get_MUA_bin(sp,WF2ephysT);
dV1 = double(dV1(1:50,:));
end