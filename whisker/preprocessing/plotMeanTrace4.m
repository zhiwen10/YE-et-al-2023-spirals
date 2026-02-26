function [h1] = plotMeanTrace4(mimgtransformed,wf_mean2)
%%
data_folder = 'E:\spiral_data_share\data';  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_25um_ssp_bfd.mat'));
atlas2 = atlas1;
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));   % atlas1
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
BW2 = BW(1:8:end,1:8:end);
%%
st2 = readtable(fullfile(data_folder,'tables',...
    'structures_ssp_bfd.csv'));                                      % a table of what all the labels mean
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
maskPath{12} = '/997/8/567/688/695/315/453/322/329/'; % SSp-bfd
maskPath{13} ='/997/8/567/688/695/315/453/322/353/'; %SSP-n
maskPath{14} ='/997/8/567/688/695/315/453/322/337/'; %SSp-ll
maskPath{15} ='/997/8/567/688/695/315/453/322/345/'; %SSp-m
maskPath{16} ='/997/8/567/688/695/315/453/322/369/'; %SSp-ul
maskPath{17} ='/997/8/567/688/695/315/453/322/361/'; %SSp-tr
maskPath{18} ='/997/8/567/688/695/315/453/322/182305689/'; %SSp-un
ssp_bfd_path = st2{st2.parent_structure_id == 329,4};
ssp_c{1} = '/997/8/567/688/695/315/453/322/329/614454341/';
ssp_c{2} = '/997/8/567/688/695/315/453/322/329/614454348/';
ssp_c{3} = '/997/8/567/688/695/315/453/322/329/614454355/';
ssp_c{4} = '/997/8/567/688/695/315/453/322/329/614454362/';
ssp_c{5} = '/997/8/567/688/695/315/453/322/329/614454369/';
ssp_c{6} = '/997/8/567/688/695/315/453/322/329/614454376/';
%%
wf1 = reshape(wf_mean2,size(wf_mean2,1)*size(wf_mean2,2),size(wf_mean2,3));
meanTrace = wf1 -mean(wf1 ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
meanTrace2 = meanTrace';
meanTrace2 = reshape(meanTrace2,size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = tracePhase';
tracePhase = reshape(tracePhase, size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));
%%
color2 = cbrewer2('seq','YlOrRd',10);
lineColor = 'k'; lineColor1 = 'w';
lineWidth1 = 1;
lineWidth2 = 1.5;
hemi = [];
scale3 = 5/8;
center = [77,104];
a = 50*2; b = 100*2; % height
c = 80*2; d = 130*2; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
ax1 = subplot(1,3,1);
mimgtransformedRGB = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = imresize(mimgtransformedRGB,[330,286]);
mimgtransformedRGB = cat(3, mimgtransformedRGB,mimgtransformedRGB,mimgtransformedRGB);
low = prctile(mimgtransformedRGB(:),2);
high = prctile(mimgtransformedRGB(:),98);
im1= imshow(mimgtransformedRGB,[low high]);
BW3 = imresize(BW,[330,286]);
set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st2,atlas2,hemi,scale3,lineColor,lineWidth1);
plotOutline(maskPath(5),st2,atlas2,[],scale3,lineColor1,lineWidth2);
for kk = 13:18
    plotOutline(maskPath(kk),st2,atlas2,[],scale3,lineColor1,lineWidth2);
end
for j = 1:39
    plotOutline(ssp_bfd_path(j),st2,atlas2,[],scale3,lineColor1,lineWidth1);
end
for j = 1:6
    plotOutline(ssp_c(j),st2,atlas2,[],scale3,[0.5,0.5,0.5],lineWidth1);
end
hold on;
th2 = 1:45:360; 
px1 = center(2);
py1 = center(1);
r = 12;
cx2 = round(r*cosd(th2)+px1);
cy2 = round(r*sind(th2)+py1);
pixel = [cy2;cx2]';
scatter(ax1,pixel(:,2)*2,pixel(:,1)*2,16,color2(3:10,:),'filled');
scatter(center(2)*2,center(1)*2,16,'w');
hold on;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
axis image; axis off;

ax2 = subplot(1,3,2);
im1= imshow(mimgtransformedRGB,[low high]);
BW3 = imresize(BW,[330,286]);
set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st2,atlas2,hemi,scale3,lineColor,lineWidth1);
plotOutline(maskPath(5),st2,atlas2,[],scale3,lineColor1,lineWidth2);
for kk = 13:18
    plotOutline(maskPath(kk),st2,atlas2,[],scale3,lineColor1,lineWidth2);
end
for j = 1:39
    plotOutline(ssp_bfd_path(j),st2,atlas2,[],scale3,lineColor1,lineWidth1);
end
for j = 1:6
    plotOutline(ssp_c(j),st2,atlas2,[],scale3,[0.5,0.5,0.5],lineWidth1);
end
hold on;
scatter(ax2,pixel(:,2)*2,pixel(:,1)*2,16,color2(3:10,:),'filled');
hold on;
scatter(center(2)*2,center(1)*2,16,'w');
axis image; axis off;
ylim([a,b]);
xlim([c,d]);

ax3 = subplot(1,3,3);
scalea = 1.5;
t1 = -2:1/35:2;
for i = 1:size(pixel,1)
    clear pks locs
    trace_raw(i,:) = wf_mean2(pixel(i,1),pixel(i,2),:);
    plot(ax3,t1(64:85),trace_raw(i,64:85)*scalea*100+scalea*(i-1),'color',color2(i+2,:),'lineWidth',1);
    % text(t1(64)-0.1,(i-1)*scalea,nameList{i},'Color', color2(i+1,:),'Interpreter','None');
    hold on; 
    [pks,locs] = findpeaks(trace_raw(i,:),'MinPeakProminence',0.001);
    scatter(t1(locs),trace_raw(i,locs)*scalea*100+scalea*(i-1),6,'k','filled');
end
ax3.FontSize = 8; 
set(ax3,'YTickLabel',[]);
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(t1(71),'--');
