function plotOverlaidAtlas(projectedAtlas2,projectedTemplate2,coords,scale)
figure;
subplot(1,2,1);
imagesc(projectedAtlas2)
axis image; axis off;
subplot(1,2,2);
imshow(projectedTemplate2,[])
axis image; axis off;
hold on;
for q = 1:numel(coords) % coords is from ctxOutlines.mat 
    cx = coords(q).x/scale;
    cy = coords(q).y/scale;    
    coordsX(q).x = cx;
    coordsX(q).y = cy;
    plot(cx,cy, 'LineWidth', 1.0, 'Color', 'r');
    hold on;
end