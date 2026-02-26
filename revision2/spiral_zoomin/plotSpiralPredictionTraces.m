function [h4bc,h4de] = plotSpiralPredictionTraces(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
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
scale3 = 5;
lineColor = 'k'; lineColor1 = 'k';
hemi = [];
%% get original and predicted dV by kernel regression.
kk = 6; % zye60
ops = get_session_info2(T,kk,data_folder);
%%
[U,V,t,mimg] = loadUVt1(ops.session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
[sp] = loadKSdir2(ops.session_root);
[syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t);
WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
dV1 = double(dV(:,~isnan(WF2ephysT1)));
V1 = double(V(1:50,~isnan(WF2ephysT1)));
[MUA_std] = get_MUA_bin(sp,WF2ephysT);
dV1 = double(dV1(1:50,:));
%%
fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
tf = load(fullfile(data_folder,'ephys','rf_tform',[fname '_tform']));
load(fullfile(data_folder,'ephys','spirals_example',[fname '_mask.mat']));
fname1 = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x'];
tf4 = load(fullfile(data_folder,'ephys','rf_tform_4x',fname1));
predict_folder = fullfile(data_folder,'ephys','dv_prediction');
load(fullfile(predict_folder,[fname '_dv_predict.mat']));
%%
epochs = [13700:13800];
frame = 45:55;
point(1,:) = [101,91]*4;
[point_t(:,1),point_t(:,2)] = transformPointsForward(tf.tform,point(:,1),point(:,2));
%%
a = 30; b = 80; % height
c = 20; d = 70; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tf.tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tf.tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tf.tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed1 = Utransformed(1:8:end,1:8:end,:);
mimgtransformed2 = mimgtransformed(1:8:end,1:8:end);
BW = logical(projectedAtlas1);
%% load dV prediction
scale = 4;
Ut = U(1:scale:end,1:scale:end,1:50);
dV1 = dV1(:,1:size(dV_predict,2));
%% angle difference for oscillation period  
[data_raw] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV1,epochs,mimgtransformed2);
[data_predict] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV_predict,epochs,mimgtransformed2);
[Axy_all] = compare_flowfield_angle(data_raw,data_predict);
%%
explained_var_all(explained_var_all(:)<0) = 0;
mean_var = squeeze(mean(explained_var_all,3));
mean_var_t = imwarp(mean_var,tf4.tform,'OutputView',imref2d(size(projectedTemplate1)));
filename = [ops.mn, '_', num2str(epochs(1))];
%%
colorn = size(MUA_std,1);
color1 = gray(colorn);
randi = randperm(colorn);    
colorAll = color1(randi,:);
%%
clear ax1 ax2 ax3 ax4
[data_raw] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV1,epochs,mimgtransformed2);
[data_predict] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV_predict,epochs,mimgtransformed2);
phase_raw = data_raw.tracePhase1;
phase_predict = data_predict.tracePhase1;
phase_raw = permute(phase_raw,[2,3,1]);
phase_predict = permute(phase_predict,[2,3,1]);
%
vx_raw = permute(data_raw.vxRaw,[2,3,1]);
vy_raw = permute(data_raw.vyRaw,[2,3,1]);
vx_predict = permute(data_predict.vxRaw,[2,3,1]);
vy_predict = permute(data_predict.vyRaw,[2,3,1]);
%% raw trace
% Ut2 = Utransformed1(:,:,1:50)./mimgtransformed2;
% Utr = reshape(Ut2,size(Ut2,1)*size(Ut2,2),size(Ut2,3));
% trace_raw1 = Utr*V(1:50,epochs);
% trace_raw = reshape(trace_raw1,size(Ut2,1),size(Ut2,2),size(trace_raw1,2));
%%
trace_raw = permute(data_raw.trace2d1,[2,3,1]);
trace_predict = permute(data_predict.trace2d1,[2,3,1]);
trace_raw2 = reshape(trace_raw,[],size(trace_raw,3));
trace_raw2(not(BW2down),:) = nan; trace_raw2 = reshape(trace_raw2, size(trace_raw));
trace_predict2 = reshape(trace_predict,[],size(trace_predict,3));
trace_predict2(not(BW2down),:) = nan; trace_predict2 = reshape(trace_predict2, size(trace_predict));
trace_raw3 = trace_raw2(a:b,c:d,:);
trace_predict3 = trace_predict2(a:b,c:d,:);
%%
row = 40, col = 40;
frame2 = 45:55;
figure;
plot([1:numel(frame2)]/35,squeeze(trace_raw3(row,col,frame2))*40,'k');
%%
hs = figure('Renderer', 'painters', 'Position', [100 100 900 300]);
ax1 = subplot(1,4,1);
scale3 = 5/8;
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
mean_var_t2 = mean_var_t(1:8:end,1:8:end);
im1 = imagesc(ax1,mean_var_t2);
set(im1, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
axis image; axis off;
cmap1 = inferno;
colormap(ax1,cmap1);
cb1 = colorbar;
cb1.Position(4) = cb1.Position(4)/2;

ax2 = subplot(1,4,2);
im2 = imagesc(ax2,mean_var_t2);
set(im2, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
axis image; axis off;
cmap1 = inferno;
colormap(ax2,cmap1);
% colorbar;
xlim([c,d]);
ylim([a,b]);

i = 3;
ax3 = subplot(1,4,3);
% frame_raw1 = squeeze(phase_raw(:,:,frame(i)));   
% im1 = imagesc(frame_raw1);
% colormap(ax3,colorcet('C06'));
% set(im1, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
% hold on;
% plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
% hold on;
% axis off; axis image;
for row = 1:4:51
    for col = 1:4:51
        plot([1:numel(frame2)]/5+col+c,squeeze(trace_raw3(row,col,frame2))*200+row+a,'k');
        hold on;
    end
end
hold on;
plot([1,10]/5+1+c, [2+row+a,2+row+a],'k');
hold on;
plot([5,5]/5+1+c, [0,0.02]*200+row+a,'k');
set(gca, 'YDir','reverse')
xlim([c,d]);
ylim([a,b]);
axis off; axis image;

ax4 = subplot(1,4,4);
% frame_predict1 = squeeze(phase_predict(:,:,frame(i)));   
% im1 = imagesc(frame_predict1);
% colormap(ax4,colorcet('C06'));
% set(im1, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
% hold on;
% plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
% axis image; axis off;
% hold on;
for row = 1:4:51
    for col = 1:4:51
        plot([1:numel(frame2)]/5+col+c,squeeze(trace_predict3(row,col,frame2))*200+row+a,'r');
        hold on;
    end
end
set(gca, 'YDir','reverse')
xlim([c,d]);
ylim([a,b]);
axis off; axis image;
%%
print(hs,'PredictionTraceExample','-dpdf', '-bestfit', '-painters');