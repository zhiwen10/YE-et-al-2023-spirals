%%
function [MatchPoint1] = UniqueMatchPoints2(points1,points2)
if size(points1,1)>size(points2)
    pointsTemp = points1;
    points1 = points2; points2 = pointsTemp;
end
% points2 with larger size
% last column has frame info
points_all = [points1;points2];
%%
[C,ia,ic] = unique(points_all(:,end));
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
MatchFrames = value_counts(value_counts(:,2)==1,1);
%%
MatchPoint1 = [];
for k = 1:numel(MatchFrames)
    pTemp1 = find(points1(:,end)==MatchFrames(k));
    pTemp2 = find(points2(:,end)==MatchFrames(k));
    if (isempty(pTemp1) & length(pTemp2)>=1)| (length(pTemp1)>=1 & isempty(pTemp2))
        points_temp = [points1(pTemp1,:);points2(pTemp2,:)];
        MatchPoint1 = [MatchPoint1;points_temp];
    end   
end