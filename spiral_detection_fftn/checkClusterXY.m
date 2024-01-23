function clustersCentroids = checkClusterXY(pwAll1,dThreshold)
clustersCentroids = [];
if not(isempty(pwAll1))
    % cluster spirals within dThreshold distance
    [~,~,clustersXY1] = clusterXYpoints(pwAll1(:,1:2),dThreshold,[],'point','merge');
    % true spirals are detected multiple times in neighboring grids
    a1 = cellfun(@(x) size(x,1)>1,clustersXY1); 
    clustersXY1 = clustersXY1(a1);
    clusterCount = numel(clustersXY1);
    for k=1:clusterCount
        clustersCentroids(k,:) = round(mean(clustersXY1{k,1},1));
    end
    % clustersCentroids(:,3) = frameN;    
    % pwAll = [pwAll; clustersCentroids];
end
end