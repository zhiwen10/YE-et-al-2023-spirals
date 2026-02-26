function [hs12g,hs12h] = plotWavePredictionExample2(T,data_folder,save_folder)
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
kk = 20; % zye67
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
epochs = [9050:9150];
frame = 44:54;
point(1,:) = [101,91]*4;
[point_t(:,1),point_t(:,2)] = transformPointsForward(tf.tform,point(:,1),point(:,2));
%%
a = 50; b = 100; % height
c = 80; d = 130; % width
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
hs12g = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
ax1 = subplot(3,4,1);
im1 = imagesc(ax1,mimgt);
set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
colormap(ax1,gray);
axis image; axis off;
mxRange = prctile(mimgt(:), 99.5);
caxis([0,mxRange])
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(point_t(:,2),point_t(:,1),6,'r','filled')
colorbar;
color1 = cbrewer2('qual','Set2',8);

ax2 = subplot(3,4,5);
im2 = imagesc(ax2,mean_var_t);
set(im2, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
axis image; axis off;
cmap1 = inferno;
colormap(ax2,cmap1);
colorbar;

point1 = round(point/scale);
mimg_scale = mimg(1:scale:end,1:scale:end);
ax3 = subplot(3,4,[2:4])
trace_predict = squeeze(data_predict.trace2d1(:,point1(1,1),point1(1,2)));
trace_raw = squeeze(data_raw.trace2d1(:,point1(1,1),point1(1,2)));
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_raw*100,'lineWidth',1,'color',color1(5,:)); % GREEN
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_predict*100,'lineWidth',1,'color',color1(8,:)); % GRAY
xlim([WF2ephysT(epochs(1)),WF2ephysT(epochs(end))]);
xlabel('Time (s)');
ylabel('df/f (%)');
hold on;
xline(WF2ephysT(epochs(frame(1))),'--');
xline(WF2ephysT(epochs(frame(end))),'--');
hold off;
%
t_sample = (epochs(end)-epochs(1));
t_length = t_sample/35;

ax4 = subplot(3,4,[6:8])
color1 = colorcet('C06','N',8);
count1 = 0;
for i = 1:numel(sp.gcluster)
    clear incl cluster_t cluster_depth
    incl = (sp.clu==sp.gcluster(i) & sp.spikeTimes>WF2ephysT(epochs(1)) ...
    & sp.spikeTimes<WF2ephysT(epochs(end)));
    cluster_t = sp.spikeTimes(incl);
    cluster_depth = sp.spikeDepths(incl);
    mean_rate = numel(cluster_t)./t_length;
    if mean_rate<20
        plot([cluster_t,cluster_t]',[cluster_depth,cluster_depth+100]','Color',colorAll(i,:),'lineWidth',1);
        hold on;
        count1 = count1+1;
    end
end
xlim([WF2ephysT(epochs(1)),WF2ephysT(epochs(end))]);

BW2down = BW2(1:8:end,1:8:end);
mean_var_t2 = mean_var_t(1:8:end,1:8:end);
indx = find(mean_var_t2<=0.4 & not(BW2down));
ax10 = subplot(3,4,10:12);
theta = angle(Axy_all);
theta(:,indx) = [];
[S s] = circ_var(theta, [], [], 2);
c1 = 1-S;
scatter(1:t_sample,c1,[],'k','filled','MarkerFaceAlpha',0.4);
xticks([0:17.5:t_sample]);
xticksN = num2cell(0:0.5:t_length);
xticksStr = cellfun(@num2str,xticksN,'UniformOutput',false);
xticklabels(xticksStr);
xlabel('Time (s)');
xlim([1,t_sample]);
ylim([0,1]);
print(hs12g, fullfile(save_folder,['FigS12g_' filename]),...
    '-dpdf', '-bestfit', '-painters');
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

v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
mean_var_t2 = mean_var_t(1:8:end,1:8:end);
BW_var = ((mean_var_t2>=0.1)& BW2down);
indx = find(not(BW_var));
theta = angle(Axy_all);
theta(:,indx) = [];
[S s] = circ_var(theta, [], [], 2);
c1 = 1-S;

clear vx_raw2 vy_raw2 vx_predict2 vy_predict2
skip = 5; zoom_scale = 2;
[vx_raw2,vy_raw2] = flow_vector_scale(vx_raw,vy_raw,skip,zoom_scale);
[vx_predict2,vy_predict2] = flow_vector_scale(vx_predict,vy_predict,skip,zoom_scale);
clear vx_raw3 vy_raw3 vx_predict3 vy_predict3
vx_predict3 = reshape(vx_predict2,[],size(vx_predict2,3));
vx_predict3(not(BW_var),:) = nan; vx_predict3 = reshape(vx_predict3, size(vx_predict2));
vy_predict3 = reshape(vy_predict2,[],size(vy_predict2,3));
vy_predict3(not(BW_var),:) = nan; vy_predict3 = reshape(vy_predict3, size(vy_predict2));
vx_raw3 = reshape(vx_raw2,[],size(vx_raw2,3));
vx_raw3(not(BW_var),:) = nan; vx_raw3 = reshape(vx_raw3, size(vx_raw2));
vy_raw3 = reshape(vy_raw2,[],size(vy_raw2,3));
vy_raw3(not(BW_var),:) = nan; vy_raw3 = reshape(vy_raw3, size(vy_raw2));
%
scale4 = 5/8;
frame_n = numel(frame);
hs12h = figure('Renderer', 'painters', 'Position', [100 100 800 400]);
for i = 1:numel(frame)   
    frame_raw = squeeze(phase_raw(:,:,frame(i)));
    
    ax1(i) = subplottight(4,frame_n,i);
    im1 = imagesc(frame_raw);
    colormap(ax1(i),colorcet('C06'));
    set(im1, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    frame_predict = squeeze(phase_predict(:,:,frame(i)));
    ax2(i) = subplottight(4,frame_n,frame_n+i);
    im2 = imagesc(frame_predict);
    colormap(ax2(i),colorcet('C06'));
    set(im2, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    ax3(i) = subplottight(4,frame_n,i+2*frame_n);
    ax3(i).Position(1) = ax3(i).Position(1)+0.005;
    ax3(i).Position(2) = ax3(i).Position(2);
    ax3(i).Position(3) = ax3(i).Position(3)-0.01;
    ax3(i).Position(4) = ax3(i).Position(4)-0.025;
    im3 = imagesc(frame_raw);
    colormap(ax3(i),colorcet('C06'));
    set(im3, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    imH1Raw4 = quiver(squeeze(vx_raw3(:,:,frame(i))),squeeze(vy_raw3(:,:,frame(i))),'k','lineWidth',0.5,'autoScale','off');   
    caxis([-pi,pi]);
    ylim([a,b]);
    xlim([c,d]);
    title(num2str(c1(frame(i))));
    
    ax4(i) = subplottight(4,frame_n,3*frame_n+i);
    ax4(i).Position(1) = ax4(i).Position(1)+0.005;
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.01;
    ax4(i).Position(4) = ax4(i).Position(4)-0.025;
    im4 = imagesc(frame_predict);
    colormap(ax4(i),colorcet('C06'));
    set(im4, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    imH1Raw4 = quiver(squeeze(vx_predict3(:,:,frame(i))),squeeze(vy_predict3(:,:,frame(i))),'k','lineWidth',0.5,'autoScale','off');  
    caxis([-pi,pi]);  
    ylim([a,b]);
    xlim([c,d]);

end
print(hs12h, fullfile(save_folder,['FigS12h_' filename '_frames']),...
    '-dpdf', '-bestfit', '-painters');