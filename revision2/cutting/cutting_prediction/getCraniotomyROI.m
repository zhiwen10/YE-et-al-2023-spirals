function [centers,radius] = getCraniotomyROI(mimgt1)
%%
centers = [];
edges = [];
figure;
ax = subplot(1,1,1);
imagesc(mimgt1);
colormap('gray');
axis image;
for k = 1:2
    roi = drawpoint(ax);
    centers(k,:) = roi.Position;
    roi = drawpoint(ax);
    edges(k,:) = roi.Position;
end
for k = 1:2
    radius(k,1) = hypot(edges(k,1)-centers(k,1),edges(k,2)-centers(k,2));
end    
centers = round(centers);
radius = round(radius);