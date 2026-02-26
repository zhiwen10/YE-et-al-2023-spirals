function h3eg = plotExampleKernelAP(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
mn = 'ZYE_0012';
td = '2020-10-16';
tdb = datestr(td,'yyyymmdd');
en = 5;
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
mimg1 = mimg/max(mimg(:));
U = U(:,:,1:50);
V = V(1:50,:);
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'tables',[fname '_tform.mat']));               % load atlas transformation matrix tform;
%%
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
%%
scale = 8;
%% right SSp index
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[Unew_SSp,Vnew_SSp] = redoSVD(UselectedSSp,V(:,1:58000));
USSp = Unew_SSp(:,1:50);
VSSp = Vnew_SSp(1:50,:);
%% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'right';
[indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[Unew_MO,Vnew_MO] = redoSVD(UselectedMO,V(:,1:58000));
UMO = Unew_MO(:,1:50);
VMO = Vnew_MO(1:50,:);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; BW_MO(indexMO) =1;
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
%% all areas
UtDown = Utransformed(1:scale:end,1:scale:end,1:50);
Utr = reshape(UtDown,size(UtDown,1)*size(UtDown,2),size(UtDown,3));
%% prepare regressor and signal for regression
regressor1 = zscore(VMO,[],2);
regressor = regressor1(:,1:58000);
signal1 = USSp*VSSp(:,1:58000);
%% regression
useGPU = 0;
if useGPU == 1
    Vregressor_GPU =gpuArray(regressor);
    signal_GPU = gpuArray(signal1);
    kk1 = gather(Vregressor_GPU'\signal_GPU');
else
    kk1 = regressor1'\signal1';
end   
k1_real = UMO*kk1;
%% prepare test regressor 
data = UselectedMO*V(:,78000:end);    
% subtract mean as we did before
data = data-mean(data,2);
regressor_test = Unew_MO' * data;   
regressor_test_z = zscore(regressor_test,[],2);
%%
projectedAtlas3 = projectedAtlas1(1:scale:end,1:scale:end);
%% prediction
predicted_signals = kk1'*regressor_test_z(1:50,:);
predicted_signal = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),size(predicted_signals,2));
predicted_signal(indexSSp,:) = predicted_signals;
rawAll = Utr*V(:,78000:end);
%% kernel project back to pixel space
kernel_full = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),size(projectedAtlas3,1)*size(projectedAtlas3,2));
for j = 1:numel(indexSSp)
     kernel_temp2 = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),1);
     kernel_temp2(indexMO) = k1_real(:,j);
     kernel_full(:,indexSSp(j)) = kernel_temp2;
end
kernel_full2 = reshape(kernel_full,[size(projectedAtlas3,1),size(projectedAtlas3,2),...
    size(projectedAtlas3,1),size(projectedAtlas3,2)]);
%%
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
point(8,:) = [110 104]; %% VISp
point(7,:) = [97 79]; %% RSP
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% 3 subareas
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
%%
pixel(1,:) = [115,105]*8; %V1
pixel(2,:) = [78,118]*8; %S1
pixel(3,:) = [89,78]*8; %RSP

% zye12 anterior cortex conterparts
pixel_a(1,:) = [64,80]*8;
pixel_a(2,:) = [54,88]*8;
pixel_a(3,:) = [72,77]*8;

% zye12 left cortex conterparts
pixel_l(1,:) = [107,41]*8;
pixel_l(2,:) = [75,30]*8;
pixel_l(3,:) = [92,61]*8;
%% color overlay a cycle
color1 = colorcet('C06', 'N', 9);
mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
h3eg = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
axx4 = subplot(2,3,1)
% im2 = imagesc(mimgtransformed2);
im2 = imagesc(template1);
colormap(gray)
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 8
TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
% overlay kernel with colormap    
maxI = max(TheColorImage(:));
[row,col] = find(TheColorImage==maxI);
imageSize = size(TheColorImage);  
colorIntensity = TheColorImage/maxI;
colorIntensity = colorIntensity.^2;
thisColor = color1(kkk,:);
thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
    thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
hold on
hb2(kkk) = imshow(thisColorImage);
hold off

set(hb2(kkk),'AlphaData',colorIntensity);
axis image; 
hold on;
scatter(axx4,col,row,12,'MarkerFaceColor','k','MarkerEdgeColor','None');
hold on;
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

axx4 = subplot(2,3,2);
im2 = imagesc(template1);
colormap(gray)
set(im2, 'AlphaData', BW_SSp, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 8;
axis image; 
hold on;
scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

Vraw = V(1:50,78000:end);
traceSSp_predict2 = reshape(predicted_signal,165,143,[]);
subplot(2,3,3);
hold off;
i = 1;
trace_raw1 = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %sensory
trace_raw3 = squeeze(Utransformed(row*scale,col*scale,:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %left hemisphere
trace_predict = squeeze(traceSSp_predict2(round(pixel(i,1)/8),round(pixel(i,2)/8),:))./mimgtransformed(pixel(i,1),pixel(i,2));  % preidicted
plot(1/35:1/35:size(trace_predict,1)/35,trace_predict+0.1,'color',color1(kkk,:));
hold on;
plot(1/35:1/35:size(trace_raw1,2)/35,trace_raw1,'color',[0 0 0 ]);
hold on;
plot(1/35:1/35:size(trace_raw3,2)/35,trace_raw3+0.2,'color',color1(kkk,:));
hold on;
xlim([50,60])
xlabel('Time (s)')
xticks([50:2:60])
xticklabels(string(num2cell(0:2:10)))
ylim([-0.05,0.25]);

% color overlay a cycle
axx4 = subplot(2,3,4);
im2 = imagesc(template1);
colormap(gray)
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 1;
TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
% overlay kernel with colormap    
maxI = max(TheColorImage(:));
[row,col] = find(TheColorImage==maxI);
imageSize = size(TheColorImage);  
colorIntensity = TheColorImage/maxI;
colorIntensity = colorIntensity.^2;
thisColor = color1(kkk,:);
thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
    thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
hold on
hb2(kkk) = imshow(thisColorImage);
hold off

set(hb2(kkk),'AlphaData',colorIntensity);
axis image; 
hold on;
scatter(axx4,col,row,12,'MarkerFaceColor','k','MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

axx4 = subplot(2,3,5)
im2 = imagesc(template1);
colormap(gray)
set(im2, 'AlphaData', BW_SSp, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 1;
axis image; 
hold on;
scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off;

Vraw = V(1:50,78000:end);
traceSSp_predict2 = reshape(predicted_signal,165,143,[]);

subplot(2,3,6);
hold off;
i = 2;
trace_raw1 = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %sensory
trace_raw3 = squeeze(Utransformed(row*scale,col*scale,:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %left hemisphere
trace_predict = squeeze(traceSSp_predict2(round(pixel(i,1)/8),round(pixel(i,2)/8),:))./mimgtransformed(pixel(i,1),pixel(i,2));  % preidicted


plot(1/35:1/35:size(trace_predict,1)/35,trace_predict+0.1,'color',color1(kkk,:));
hold on;
plot(1/35:1/35:size(trace_raw1,2)/35,trace_raw1,'color',[0 0 0 ]);
hold on;
plot(1/35:1/35:size(trace_raw3,2)/35,trace_raw3+0.2,'color',color1(kkk,:));

xlim([50,60])
xlabel('Time (s)')
xticks([50:2:60])
xticklabels(string(num2cell(0:2:10)))
ylim([-0.05, 0.25]);
%%
print(h3eg, fullfile(save_folder,'Fig3eg_example_traces_MO2SSp'),...
    '-dpdf', '-bestfit', '-painters');