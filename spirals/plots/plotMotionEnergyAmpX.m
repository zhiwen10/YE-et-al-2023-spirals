function hs8gh = plotMotionEnergyAmpX(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
scale = 1;
%%
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
areaPath{1} = '/997/8/567/688/695/315/453/322/';                           % SSp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
scale = 1;
[indexright,UselectedRight] = select_area(sensoryArea,spath,...
    st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_right = BW_empty; 
BW_right(indexright) =1;
[row,col] = find(BW_right);
brain_index = [col,row];
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
edges = [0:100000:1000000];
%%
frame_all = 0;
hist_bin = 40;
max_intensity = nan(14,10);
count1 = 1;
R1 = []; p1 = [];
for kk = [1:6,9:15]
    %%
    clear indx2 spiralsT lia locb
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,'spirals','spirals_index',[fname '_amp']));
    load(fullfile(data_folder,'spirals','spirals_index',[fname '_motion_energy']));
    tsize = min(size(traceAmp_mean,1), size(image_energy2,1));
    traceAmp_mean = traceAmp_mean(1:tsize,1); 
    image_energy2 = image_energy2(1:tsize,1);
    data1 = [traceAmp_mean image_energy2];
    indx1 = not(isnan(data1(:,2)));
    data1 = data1(indx1,:);
    %%
    [R,p] = corrcoef(data1(:,1),data1(:,2));
    R1(count1,1) = R(2,1);
    p1(count1,1) = p(2,1);
    %%
    count1 = count1+1;
end
%%
indx = randperm(size(data1,1),10000);
data2 = data1(indx,:);
mdl = fitlm(data2(:,1),data2(:,2));
x1 = 0.002:0.001:0.015;
y1 = mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*x1;
hs8gh = figure; 
subplot(1,2,1);
scatter_kde(data2(:,1),data2(:,2),'filled','MarkerSize', 3);
colormap(hot);
hold on;
plot(x1,y1,'k')
xlabel('2-8Hz amplitude');
ylabel('Motion energy');
text(0.005,1400000,['r = ' num2str(round(R1(end)*10)/10)]);
subplot(1,2,2);
hold off;
mean_r = mean(R1);
sem_r = std(R1)/sqrt(13);
rnd1 = normrnd(0,0.1,[1,13]);
xa = zeros(1,13)+rnd1;
bar(0,mean_r);
hold on;
scatter(xa,R1);
hold on;
errorbar(0,mean_r,sem_r,'k');
xlim([-1,1]);
ylim([-1,0.2]);
%%
print(hs8gh, fullfile(save_folder,'Figs8gh_motion_amp_correlation.pdf'),...
    '-dpdf', '-bestfit', '-painters');