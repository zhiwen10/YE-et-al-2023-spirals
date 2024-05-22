
function probeSites = makeProbePlot(ax, cm, ss)
% ax is figure axis
% cm is chanMap struture loaded from .mat file
% ss controls the size of the sites plotted
hold(ax, 'on');
nSites = numel(cm.xcoords);

probeSites = [];
% 4 conners of the site square.
sq = ss*([0 0; 0 1; 1 1; 1 0]-[0.5 0.5]);
% add sites coordinates to the square, and draw square polygons
for q = 1:nSites
    probeSites(q) = fill(ax, ...
        sq(:,1)+cm.xcoords(q), ...
        sq(:,2)+cm.ycoords(q), 'b', 'EdgeAlpha',0);    
end


