function overlayOutlines(coords,scale,color)
if nargin>2
    color1 = color;
else
    color1 = 'k';
end
for q = 1:numel(coords) % coords is from ctxOutlines.mat 
    cx = coords(q).x/scale;
    cy = coords(q).y/scale;    
    coordsX(q).x = cx;
    coordsX(q).y = cy;
    plot(cx,cy, 'LineWidth', 0.5, 'Color', color1);
    hold on;
end
end