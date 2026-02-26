function [U, V, t, mimg] = loadUVt1(expRoot)

U = readNPY(fullfile(expRoot,'svdSpatialComponents.npy'));
mimg = readNPY(fullfile(expRoot, 'meanImage.npy'));
V = readNPY(fullfile(expRoot,'svdTemporalComponents_corr.npy'));
t = readNPY(fullfile(expRoot,'svdTemporalComponents_corr.timestamps.npy'));

if length(t) > size(V,2)
    t = t(1:size(V,2));
elseif length(t)< size(V,2)
    V = V(:,1:numel(t));
end