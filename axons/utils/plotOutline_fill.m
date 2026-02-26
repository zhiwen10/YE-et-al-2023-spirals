function plotOutline_fill(areaPath,st,section,hemisphere,scale,fillcolor, faceAlpha)
% cortex path
% areaPath = '/997/8/567/623/477/';
if nargin == 6
    faceAlpha = 1;
end
    
indx = [];
for r = 1:numel(areaPath)
    areaPath1 = areaPath{r};
    indx1 = find(cellfun(@(x)contains(x,areaPath1),st.structure_id_path));
    indx = [indx; indx1];
end

idAll = st.id(indx);
section = imresize(section,scale);
indx2 = ismember(double(section(:)),double(idAll));
area = zeros(size(section));
area(indx2) = 1;
if hemisphere == -1 % plot left
    area(:,size(section,2)/2+1:size(section,2)) = 0;
elseif hemisphere == 1 % plot right
    area(:,1:size(section,2)/2) = 0;
end
area1 = logical(area);
c1 = contourc(double(area1>0), [0.5 0.5]);
coordsReg1 = makeSmoothCoords(c1);
for cidx1 = 1:numel(coordsReg1)
    h = fill(gca, coordsReg1(cidx1).x,coordsReg1(cidx1).y, fillcolor, 'EdgeColor', 'k','FaceAlpha',faceAlpha,'lineWidth',1); 
    % h = plot(gca, coordsReg1(cidx1).x*scale,coordsReg1(cidx1).y*scale, 'k'); 
    h.LineWidth = 1.0;
    hold on;
end
end