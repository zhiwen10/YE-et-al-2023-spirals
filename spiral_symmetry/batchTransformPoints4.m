%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
mn1{1} = 'AB_0004'; td1{1} = '2021-03-30'; en1{1} = 1;
mn1{2} = 'ZYE_0052'; td1{2} = '2021-12-18'; en1{2} = 2;
mn1{3} = 'ZYE_0056'; td1{3} = '2022-01-10'; en1{3} = 1;
mn1{4} = 'ZYE_0058'; td1{4} = '2022-03-09'; en1{4} = 1;
mn1{5} = 'ZYE_0012'; td1{5} = '2020-10-16'; en1{5} = 5;
%%
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
spiralAll = [];
colorAll = [];
%
for k  = 1:5
    %%
    clearvars -except mn1 td1 en1 folder spiralAll colorAll k 
    %%
    mn = mn1{k}; td = td1{k}; en = en1{k};
    formatOut = 'yyyymmdd';
    tdb = datestr(td,formatOut);
    fname = [mn '_' tdb '_' num2str(en)];
    dfolder = fullfile(folder,mn,td,num2str(en));
    load(fullfile(dfolder,[fname '_filtered_spirals.mat']));
    dfolder = fullfile(folder,mn1{k},td1{k},num2str(en1{k}));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    
    filteredSpirals2(filteredSpirals2(:,4) == -1,4) = 0;
    spirals2 = filteredSpirals2;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,spirals2(:,1),spirals2(:,2));
    color1 = zeros(size(spiralsT ,1),3);
    color1(:,1) = spirals2(:,end-1);
    color1(:,2) = double(not(spirals2(:,end-1)));
    spirals2(:,5) = spirals2(:,5)+k*10^8;
    spirals2(:,1:2) = spiralsT;    
    spiralAll = [spiralAll;spirals2];
    colorAll = [colorAll;color1];
end
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
%%
h1 = figure;
scatter(spiralAll(1:50000,1),spiralAll(1:50000,2),3,colorAll(1:50000,:),'filled');
set(gca, 'ydir','reverse')
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
%%
figure;
h = scatter_kde(spiralAll(:,1),spiralAll(:,2),'filled','MarkerSize',3);
set(gca, 'ydir','reverse')
axis image
legend off;
colormap(hot)
hold on;
scale2=1;
overlayOutlines(coords,scale2);
axis off;
%%
roiM1L = drawpolygon; roiSL = drawpolygon; roiM1R = drawpolygon; 
roiSR = drawpolygon; roiM2L = drawpolygon; roiM2R = drawpolygon; 
%%
load('roiSelection.mat')
%%
for k  = 1:5
    %%
    clear spiralsT spirals2
    %%
    mn = mn1{k}; td = td1{k}; en = en1{k};
    formatOut = 'yyyymmdd';
    tdb = datestr(td,formatOut);
    fname = [mn '_' tdb '_' num2str(en)];
    dfolder = fullfile(folder,mn,td,num2str(en));
    load(fullfile(dfolder,[fname '_filtered_spirals.mat']));
    dfolder = fullfile(folder,mn1{k},td1{k},num2str(en1{k}));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    
    filteredSpirals2(filteredSpirals2(:,4) == -1,4) = 0;
    spirals2 = filteredSpirals2;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,spirals2(:,1),spirals2(:,2));
    color1 = zeros(size(spiralsT ,1),3);
    color1(:,1) = spirals2(:,end-1);
    color1(:,2) = double(not(spirals2(:,end-1)));
    spirals2(:,5) = spirals2(:,5)+k*10^8;
    spirals2(:,1:2) = spiralsT;    
    ratio = get_ratio(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spirals2);
    ratio_all(:,k) = ratio;
end
%%
area_name = {'M1-S1-left','M1-S1-right','S1-left-right','M1-left-right','M1-M2-left','M1-M2-right'};
T = array2table(ratio_all,'RowNames',area_name,'VariableNames',mn1);
T.mean = mean(ratio_all,2);
T.std = std(ratio_all,[],2);
r = normrnd(0,0.05,[1,5]);

h1 = figure; 
bar(T.mean,'faceColor',[0.5, 0.5, 0.5])
for i = 1:6
    hold on;
    scatter(ones(5,1)'*i+r,T{i,1:5},'filled');
end
xticks([1 2 3 4 5 6]);
xticklabels(area_name);
ylabel('matching ratio');
%%
print(h1, 'spiral symmetry ratio', '-dpdf', '-bestfit', '-painters');
%%
figure; 
bar(T.mean(2:3),'faceColor',[0.5, 0.5, 0.5]);
for i = 2:3
    hold on;
    scatter(ones(5,1)'*(i-1)+r,T{i,1:5},'filled')
end
xticks([1 2]);
xticklabels(area_name(2:3));
ylabel('matching ratio');

%%
function ratio = get_ratio(roiM1L,roiSL,roiM1R,roiSR,roiM2L,roiM2R,spiralAll)
%%
tf1 = inROI(roiM1L,spiralAll(:,1),spiralAll(:,2)); pointsM1L = spiralAll(tf1,:); 
tf2 = inROI(roiSL,spiralAll(:,1),spiralAll(:,2)); pointsSL = spiralAll(tf2,:);
tf3 = inROI(roiM1R,spiralAll(:,1),spiralAll(:,2)); pointsM1R = spiralAll(tf3,:);
tf4 = inROI(roiSR,spiralAll(:,1),spiralAll(:,2)); pointsSR = spiralAll(tf4,:);
tf5 = inROI(roiM2L,spiralAll(:,1),spiralAll(:,2)); pointsM2L = spiralAll(tf5,:);
tf6 = inROI(roiM2R,spiralAll(:,1),spiralAll(:,2)); pointsM2R = spiralAll(tf6,:);
%%
[MatchPoint1,MatchPoint2,DirDiffa] = UniqueMatchPoints(pointsM1L,pointsSL);
Indxa = find(DirDiffa~=0); color2 = zeros(size(MatchPoint1,1),1); color2(Indxa,1) = 1;
colora = zeros(size(MatchPoint1,1),1); colora(Indxa,1) = 1;
% find matched points
frameSML1 = MatchPoint1(DirDiffa==0,:);
frameSML2 = MatchPoint2(DirDiffa==0,:);
[MatchPoint3,MatchPoint4,DirDiffb] = UniqueMatchPoints(pointsM1R,pointsSR);
Indxb = find(DirDiffb~=0);
colorb = zeros(size(MatchPoint3,1),1); colorb(Indxb,1) = 1;
frameSMR1 = MatchPoint3(DirDiffb==0,:);
frameSMR2 = MatchPoint4(DirDiffb==0,:);
%%
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
[MatchPoint7,MatchPoint8,DirDiffd] = UniqueMatchPoints(pointsM1L,pointsM1R);
Indxd = find(DirDiffd~=0);
colord = zeros(size(MatchPoint7,1),1); colord(Indxd,1) = 1;

[MatchPoint9,MatchPoint10,DirDiffe] = UniqueMatchPoints(pointsM2L,pointsM2R);
Indxe = find(DirDiffe~=0);
colore = zeros(size(MatchPoint9,1),1); colore(Indxe,1) = 1;
%%
[MatchPoint11,MatchPoint12,DirDifff] = UniqueMatchPoints(pointsM1L,pointsM2L);
Indxf = find(DirDifff~=0); 
colorf = zeros(size(MatchPoint11,1),1); colorf(Indxf,1) = 1;
% find matched points
frameM12L = MatchPoint11(DirDifff==0,:);

[MatchPoint13,MatchPoint14,DirDiffg] = UniqueMatchPoints(pointsM1R,pointsM2R);
Indxg = find(DirDiffg~=0);
colorg = zeros(size(MatchPoint13,1),1); colorg(Indxg,1) = 1;
%%
num1  =5; num2  =2;
color_new = cbrewer2('qual','Set2',8);
colora1 = repmat(color_new(num1,:),size(colora,1),1);
colora1(logical(colora),:) = repmat(color_new(num2,:),sum(colora),1);

colorb1 = repmat(color_new(num1,:),size(colorb,1),1);
colorb1(logical(colorb),:) = repmat(color_new(num2,:),sum(colorb),1);

colorc1 = repmat(color_new(num1,:),size(colorc,1),1);
colorc1(logical(colorc),:) = repmat(color_new(num2,:),sum(colorc),1);

colord1 = repmat(color_new(num1,:),size(colord,1),1);
colord1(logical(colord),:) = repmat(color_new(num2,:),sum(colord),1);

colore1 = repmat(color_new(num1,:),size(colore,1),1);
colore1(logical(colore),:) = repmat(color_new(num2,:),sum(colore),1);

colorf1 = repmat(color_new(num1,:),size(colorf,1),1);
colorf1(logical(colorf),:) = repmat(color_new(num2,:),sum(colorf),1);

colorg1 = repmat(color_new(num1,:),size(colorg,1),1);
colorg1(logical(colorg),:) = repmat(color_new(num2,:),sum(colorg),1);
%%
LM1Smatch = numel(find(colora==0)); LM1Sno = numel(find(colora==1)); LM1Ssum = numel(colora);
RM1Smatch = numel(find(colorb==0)); RM1Sno = numel(find(colorb==1)); RM1Ssum = numel(colorb);
LRSmatch = numel(find(colorc==0)); LRSno = numel(find(colorc==1)); LRSsum = numel(colorc);
LRM1match = numel(find(colord==0)); LRM1no = numel(find(colord==1)); LRM1sum = numel(colord);
LM12match = numel(find(colorf==0)); LM12no = numel(find(colorf==1)); LM12sum = numel(colorf);
RM12match = numel(find(colorg==0)); RM12no = numel(find(colorg==1)); RM12sum = numel(colorg);

name = {'M1-S1-left','M1-S1-right','S1-left-right','M1-left-right','M1-M2-left','M1-M2-right'};
LM1Ssum = (LM1Smatch/LM1Ssum);
RM1Ssum = (RM1Smatch/RM1Ssum);
LRSsum = (LRSmatch/LRSsum);
LRM1sum = (LRM1match/LRM1sum);
LM12sum = (LM12match/LM12sum);
RM12sum = (RM12match/RM12sum);
ratio = [LM1Ssum;RM1Ssum;LRSsum;LRM1sum;LM12sum;RM12sum];
end
