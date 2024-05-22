function colorSites(probeSites, data, colMap, cax)

colors = zeros(numel(data), 3); 
colVals = linspace(cax(1),cax(2),size(colMap,1));
data(data<=cax(1)) = cax(1);
data(data>=cax(2)) = cax(2);
for c = 1:3
    colors(:,c) = interp1(colVals, colMap(:,c), data);
end
for n = 1:numel(probeSites)
    set(probeSites(n), 'FaceColor', colors(n,:))
end

