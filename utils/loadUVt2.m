function [U, V, t, mimg] = loadUVt2(expRoot)

load(fullfile(expRoot,'svdSpatialComponents.mat')); % U
load(fullfile(expRoot, 'meanImage.mat')); %mimg
load(fullfile(expRoot,'svdTemporalComponents_corr.mat')); % V
load(fullfile(expRoot,'svdTemporalComponents_corr_timestamps.mat')); %t

if length(t) > size(V,2)
    t = t(1:size(V,2));
elseif length(t)< size(V,2)
    V = V(:,1:numel(t));
end