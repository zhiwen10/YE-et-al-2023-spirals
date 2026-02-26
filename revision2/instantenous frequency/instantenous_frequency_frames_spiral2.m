githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
pixSize = 3.45/1000/0.6*3;  % mm / pix for spiral radius
%%
scale1 = 16;
th2 = 1:5:360; 
frt_all = [];
%%
for kk = 1:size(T,1)
    %%
    clear nframe frt filteredSpirals2 meanTrace
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    %%
    filteredSpirals2 = []; 
    for m = 1:size(archiveCell,1)
        temp1 = archiveCell{m};
        slength = size(temp1,1);
        temp1(:,6) = slength;
        filteredSpirals2 = [filteredSpirals2;temp1];
    end
    filteredSpirals2(filteredSpirals2(:,4)== 0,4) = -1;
    filteredSpirals2 = filteredSpirals2(filteredSpirals2(:,5)<60000,:);
    [filteredSpirals2(:,1),filteredSpirals2(:,2)] = transformPointsForward(...
    tform,filteredSpirals2(:,1),filteredSpirals2(:,2));  
    filteredSpirals2(:,1:2) = round(filteredSpirals2(:,1:2)); 
    indx1 = (filteredSpirals2(:,1)<950 & filteredSpirals2(:,1)>850 ...
        & filteredSpirals2(:,2)<700 & filteredSpirals2(:,2)>500);
    filteredSpirals2 = filteredSpirals2(indx1,:);
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)]; 
    %%
    Ut = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Ut2 = Ut(1:16:end,1:16:end,1:50);
    Utr = reshape(Ut2,size(Ut2,1)*size(Ut2,2),size(Ut2,3));
    %%
    freq = [2,8];
    Fs = 35;
    trace = Utr*dV(1:50,1:60000);
    meanTrace = trace -mean(trace,2);
    meanTrace = double(meanTrace)';             
%     [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
%     meanTrace = filtfilt(f1,f2,meanTrace);
%%
    traceHilbert =hilbert(meanTrace);
    tracePhase = angle(traceHilbert);
    tracePhase = reshape(tracePhase,[],size(Ut2,1),size(Ut2,2));
    tracePhase = permute(tracePhase,[2,3,1]);
    %%    
    for m = 1:size(filteredSpirals2,1) 
        %%
        clear tracePhaseA tracePhaseB ft cx3 cy3
        cx2 = []; cy2 = [];
        spiral_temp = filteredSpirals2(m,:);
        frame_n = spiral_temp(1,5);
        px1 = spiral_temp(1,1)/scale1;
        py1 = spiral_temp(1,2)/scale1;
        r = spiral_temp(1,3)*pixSize/(0.01*scale1);
        for r1 = 1:r
            cx2 = [cx2,round(r1*cosd(th2)+px1)];
            cy2 = [cy2,round(r1*sind(th2)+py1)];
        end
        %%
%         figure;        
%         ax1 = subplot(1,1,1);
%         im_phase = imagesc(tracePhase(:,:,frame_n));
%         colormap(ax1,colorcet('C06'))
%         axis image; axis off;
%         for r1 = 1:r
%             cx3 = r1*cosd(th2)+px1;
%             cy3 = r1*sind(th2)+py1;
%             hold on;
%             if spiral_temp(1,4) == 1                                               % counterclockwise, then color white
%                 color1 = 'w';
%             else
%                 color1 = 'k';                                                      % clockwise, then color black
%             end
%             plot([cx3 cx3(1)],[cy3 cy3(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
%             hold on;
%             scatter(px1,py1,8,color1,'filled');
%         end
        %%
        for k = 1:length(cx2)
            tracePhaseA(k,1) = tracePhase(cy2(k),cx2(k),frame_n);
            tracePhaseB(k,1) = tracePhase(cy2(k),cx2(k),frame_n+1);
        end
        ft = angle(mean(exp(i*(tracePhaseB-tracePhaseA))));
        ft = abs(ft);
        frt(m,1) = ft*35/(2*pi);
        frt(m,2) = spiral_temp(1,6);
    end 
    frt_all{kk} = frt;
end
%%
save('instantenous_frequency_across_sessions_dV2.mat','frt_all');

%%
frt_all1 = vertcat(frt_all{:});
frt_all1(:,2) = frt_all1(:,2)/35;
h1a = figure('Renderer', 'painters', 'Position', [100 100 300 100]);
ax1 = subplot(1,1,1);
scatter_kde(frt_all1(:,1), frt_all1(:,2), 'filled', 'MarkerSize', 32);
xlim([0,15]);
box on;
% colormap(ax1,flipud(gray));
% Add Color bar
cb = colorbar();
cb.Label.String = 'Density';
print(h1a, 'Instantenous frequency_scatter.pdf','-dpdf', '-bestfit', '-painters');
%% 
clear ratio1 N1
edges = 1:50;
for n = 1:15
    frt1 = frt_all{n};
    counts1 = histcounts(frt1(:,2),edges);
    N1(:,n) = counts1;
    N_ratio_all1(:,n) = N1(:,n)./size(frt1,1);
end
mean_N = mean(N_ratio_all1,2);
std_N = std(N_ratio_all1,[],2);
h1d = figure('Renderer', 'painters', 'Position', [100 100 600 200]);
subplot(1,3,1);
x = [1:49]/35;
bar(x,mean_N,1,'w');
hold on;
errorbar(x,mean_N, std_N,'r',"CapSize",6,"LineStyle","none");
xlim([0,15/35]);
ylim([0,1]);
xticks([0, 0.2, 0.4, 0.6,0.8 1.0])
t1 = [0, 0.2, 0.4, 0.6,0.8 1.0];
t2 = cellstr(string(t1));
xticklabels(t2);
xlabel('Spiral duration (ms)');
ylabel('Spiral ratio');

% frequency
subplot(1,3,2);
clear ratio1 N1
edges = 0:0.5:15;
for n = 1:15
    frt1 = frt_all{n};
    counts1 = histcounts(frt1(:,1),edges);
    N1(:,n) = counts1;
    ratio1(:,n) = N1(:,n)./size(frt1,1);
end
ratio_mean = mean(ratio1,2);
ratio_sem = std(ratio1,[],2)./sqrt(15);
% s1 = rectangle('Position',[2 0 6 0.25],'FaceColor',[0.5 .5 .5]);
bar(edges(1:end-1),ratio_mean, 1,'w');
hold on;
s2 = errorbar(edges(1:end-1),ratio_mean,ratio_sem,'r',"CapSize",4,"LineStyle","none");
hold on;
[ma,mindex] = max(ratio_mean);
s3 = xline(edges(mindex),'--');
text(edges(mindex),0.08,['maxFreq = ' num2str(edges(mindex)) 'Hz']);
yticks([0:0.05:0.25]);
ylim([0,0.25]);
yticklabels({[0:0.05:0.25]});
xlim([0,15]);
xlabel('Instantenous frequency (Hz)');
ylabel('Proportion of spirals');
set(gca, 'box', 'off');

subplot(1,3,3);
for n = 1:15
    clear frt
    frt = frt_all{n};
    for mm = 1:14
        clear frt_temp
        frt_temp = frt(frt(:,2)==mm,1);
        frt_mean(n,mm) = mean(frt_temp,"omitnan");
        frt_sem(n,mm) = std(frt_temp,[],1,"omitnan")./sqrt(length(frt_temp));
    end
end
frt_mean2 = mean(frt_mean,1,"omitnan");
count1 = sum(not(isnan(frt_mean)),1);
frt_std2 = std(frt_mean,[],1,"omitnan");
for m = 1:14
    frt_sem2(m,1) = frt_std2(m)./sqrt(count1(m));
end
errorbar([1:10]/35,frt_mean2(1:10),frt_sem2(1:10),'k');
xlim([0,0.3]);
ylim([4,5]);
% ylim([3,8]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');
%%
print(h1d, 'Instantenous frequency_filtered2.pdf','-dpdf', '-bestfit', '-painters');