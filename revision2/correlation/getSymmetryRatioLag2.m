function ratio = getSymmetryRatioLag2(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spiralAll,lag)
%%
tf2 = inROI(roiSL,spiralAll(:,1),spiralAll(:,2)); pointsSL = spiralAll(tf2,:);
pointsM1L = pointsSL;
pointsM1L(:,end) = pointsSL(:,end)+lag;
%%
% left M1 and S1
[MatchPoint1,MatchPoint2,DirDiffa] = UniqueMatchPoints(pointsM1L,pointsSL);
Indxa = find(DirDiffa~=0); color2 = zeros(size(MatchPoint1,1),1); color2(Indxa,1) = 1;
colora = zeros(size(MatchPoint1,1),1); colora(Indxa,1) = 1;
% find matched points
frameSML1 = MatchPoint1(DirDiffa==0,:);
frameSML2 = MatchPoint2(DirDiffa==0,:);
%%
num1  =5; num2  =2;
color_new = cbrewer2('qual','Set2',8);
colora1 = repmat(color_new(num1,:),size(colora,1),1);
colora1(logical(colora),:) = repmat(color_new(num2,:),sum(colora),1); % left M1 and S1
%%
LM1Smatch = numel(find(colora==0)); LM1Sno = numel(find(colora==1)); LM1Ssum = numel(colora);
name = {'S1-left-lag'};
LM1Ssum = (LM1Smatch/LM1Ssum);
ratio = [LM1Ssum];
end
