function ratio = getSymmetryRatioLag(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spiralAll,lag)
%%
tf1 = inROI(roiM1L,spiralAll(:,1),spiralAll(:,2)); pointsM1L = spiralAll(tf1,:); 
tf2 = inROI(roiSL,spiralAll(:,1),spiralAll(:,2)); pointsSL = spiralAll(tf2,:);
tf3 = inROI(roiM1R,spiralAll(:,1),spiralAll(:,2)); pointsM1R = spiralAll(tf3,:);
tf4 = inROI(roiSR,spiralAll(:,1),spiralAll(:,2)); pointsSR = spiralAll(tf4,:);
tf5 = inROI(roiM2L,spiralAll(:,1),spiralAll(:,2)); pointsM2L = spiralAll(tf5,:);
tf6 = inROI(roiM2R,spiralAll(:,1),spiralAll(:,2)); pointsM2R = spiralAll(tf6,:);
%%
pointsSL(:,end) = pointsSL(:,end)+lag;
%%
% left M1 and S1
[MatchPoint1,MatchPoint2,DirDiffa] = UniqueMatchPoints(pointsM1L,pointsSL);
Indxa = find(DirDiffa~=0); color2 = zeros(size(MatchPoint1,1),1); color2(Indxa,1) = 1;
colora = zeros(size(MatchPoint1,1),1); colora(Indxa,1) = 1;
% find matched points
frameSML1 = MatchPoint1(DirDiffa==0,:);
frameSML2 = MatchPoint2(DirDiffa==0,:);

% right M1 and S1
[MatchPoint3,MatchPoint4,DirDiffb] = UniqueMatchPoints(pointsM1R,pointsSR);
Indxb = find(DirDiffb~=0);
colorb = zeros(size(MatchPoint3,1),1); colorb(Indxb,1) = 1;
frameSMR1 = MatchPoint3(DirDiffb==0,:);
frameSMR2 = MatchPoint4(DirDiffb==0,:);
%%
% S1 left and right
fullXSize = 1200;
symetryCriteria = 100;
[MatchPoint5,MatchPoint6,DirDiffc] = UniqueMatchPoints(pointsSL,pointsSR);
xcheck1 = MatchPoint5(:,1)+MatchPoint6(:,1); ind1 = find(xcheck1>=(fullXSize-symetryCriteria) & xcheck1<=(fullXSize+symetryCriteria));
ycheck1 = MatchPoint5(:,2)-MatchPoint6(:,2); ind2 = find(ycheck1>=-symetryCriteria & ycheck1<=symetryCriteria);
ind = intersect(ind1,ind2);
MatchPoint5 = MatchPoint5(ind,:); MatchPoint6 = MatchPoint6(ind,:); DirDiffc = DirDiffc(ind,:);
Indxc = find(DirDiffc~=0);
colorc = zeros(size(MatchPoint5,1),1); colorc(Indxc,1) = 1;
% get frame number of un-matching left-right S1
frameS1 = MatchPoint5(Indxc,end);
%%
% M1 left and right
[MatchPoint7,MatchPoint8,DirDiffd] = UniqueMatchPoints(pointsM1L,pointsM1R);
Indxd = find(DirDiffd~=0);
colord = zeros(size(MatchPoint7,1),1); colord(Indxd,1) = 1;

% M2 left and right
[MatchPoint9,MatchPoint10,DirDiffe] = UniqueMatchPoints(pointsM2L,pointsM2R);
Indxe = find(DirDiffe~=0);
colore = zeros(size(MatchPoint9,1),1); colore(Indxe,1) = 1;
%%
% left M1 and M2
[MatchPoint11,MatchPoint12,DirDifff] = UniqueMatchPoints(pointsM1L,pointsM2L);
Indxf = find(DirDifff~=0); 
colorf = zeros(size(MatchPoint11,1),1); colorf(Indxf,1) = 1;
% find matched points
frameM12L = MatchPoint11(DirDifff==0,:);

% right M1 and M2
[MatchPoint13,MatchPoint14,DirDiffg] = UniqueMatchPoints(pointsM1R,pointsM2R);
Indxg = find(DirDiffg~=0);
colorg = zeros(size(MatchPoint13,1),1); colorg(Indxg,1) = 1;
%%
num1  =5; num2  =2;
color_new = cbrewer2('qual','Set2',8);
colora1 = repmat(color_new(num1,:),size(colora,1),1);
colora1(logical(colora),:) = repmat(color_new(num2,:),sum(colora),1); % left M1 and S1

colorb1 = repmat(color_new(num1,:),size(colorb,1),1);
colorb1(logical(colorb),:) = repmat(color_new(num2,:),sum(colorb),1); % right M1 and S1

colorc1 = repmat(color_new(num1,:),size(colorc,1),1);
colorc1(logical(colorc),:) = repmat(color_new(num2,:),sum(colorc),1); % S1 left and right

colord1 = repmat(color_new(num1,:),size(colord,1),1);
colord1(logical(colord),:) = repmat(color_new(num2,:),sum(colord),1); % M1 left and right

colore1 = repmat(color_new(num1,:),size(colore,1),1);
colore1(logical(colore),:) = repmat(color_new(num2,:),sum(colore),1); % M2 left and right

colorf1 = repmat(color_new(num1,:),size(colorf,1),1);
colorf1(logical(colorf),:) = repmat(color_new(num2,:),sum(colorf),1); % left M1 and M2

colorg1 = repmat(color_new(num1,:),size(colorg,1),1);
colorg1(logical(colorg),:) = repmat(color_new(num2,:),sum(colorg),1); % right M1 and M2
%%
LM1Smatch = numel(find(colora==0)); LM1Sno = numel(find(colora==1)); LM1Ssum = numel(colora);
RM1Smatch = numel(find(colorb==0)); RM1Sno = numel(find(colorb==1)); RM1Ssum = numel(colorb);
LRSmatch = numel(find(colorc==0)); LRSno = numel(find(colorc==1)); LRSsum = numel(colorc);
LRM1match = numel(find(colord==0)); LRM1no = numel(find(colord==1)); LRM1sum = numel(colord);
LRM2match = numel(find(colore==0)); LRM2no = numel(find(colore==1)); LRM2sum = numel(colore);
LM12match = numel(find(colorf==0)); LM12no = numel(find(colorf==1)); LM12sum = numel(colorf);
RM12match = numel(find(colorg==0)); RM12no = numel(find(colorg==1)); RM12sum = numel(colorg);

name = {'M1-S1-left','M1-S1-right','S1-left-right','M1-left-right','M2-left-right','M1-M2-left','M1-M2-right'};
LM1Ssum = (LM1Smatch/LM1Ssum);
RM1Ssum = (RM1Smatch/RM1Ssum);
LRSsum = (LRSmatch/LRSsum);
LRM1sum = (LRM1match/LRM1sum);
LRM2sum = (LRM2match/LRM2sum);
LM12sum = (LM12match/LM12sum);
RM12sum = (RM12match/RM12sum);
ratio = [LM1Ssum;RM1Ssum;LRSsum;LRM1sum;LRM2sum;LM12sum;RM12sum];
end
