function [U,V,t,mimg] = get_wf_svd1(serverRoot)
%plot PCA components and traces for blue and purple channels
corrPath = fullfile(serverRoot, 'corr', 'svdTemporalComponents_corr.npy');
if ~exist(corrPath, 'file')
    %%
    colors = {'blue', 'violet'};
    computeWidefieldTimestamps(serverRoot, colors); % preprocess video
    nSV = 200;
    [U, V, t, mimg] = hemoCorrect(serverRoot, nSV); % process hemodynamic correction
else
    nSV = 200;
    [U, V, t, mimg] = loadUVt(serverRoot, nSV);
end
if length(t) > size(V,2)
  %t = t(1:end-1);
    t = t(1:size(V,2));
elseif length(t)< size(V,2)
    V = V(:,1:numel(t));
end
%%
% pixelCorrelationViewerSVD(U,V)
% dV = [zeros(size(V,1),1) diff(V,[],2)];
% ddV = [zeros(size(dV,1),1) diff(dV,[],2)];