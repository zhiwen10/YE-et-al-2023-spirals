%%
function [MatchPoint1,MatchPoint2,DirDiff] = UniqueMatchPoints2(points1,points2)
if size(points1,1)>size(points2)
    pointsTemp = points1;
    points1 = points2; points2 = pointsTemp;
end
% points2 with larger size
% last column has frame info
[C,ia,ic] = unique(points2(:,end));
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
singlePoints2Frames = value_counts(value_counts(:,2)==1,1);
[fa,fb] = ismember(singlePoints2Frames,points1(:,end));
MatchFrames = singlePoints2Frames(fa);
%%
kk=1;
for k = 1:numel(MatchFrames)
    pTemp1 = find(points1(:,end)==MatchFrames(k));
    pTemp2 = find(points2(:,end)==MatchFrames(k));
    if length(pTemp1) == 1 && length(pTemp2) == 1
        MatchPoint1(kk,:) = points1(pTemp1,:);
        MatchPoint2(kk,:) = points2(pTemp2,:);
        kk = kk+1;
    end   
end
DirPoints1 = MatchPoint1(:,end-1); DirPoints2 = MatchPoint2(:,end-1);
DirPoints1Invert =  -DirPoints1;
DirDiff = DirPoints1Invert-DirPoints2;