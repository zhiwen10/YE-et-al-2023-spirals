%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%% only select cortex in the atlas
spath = string(st.structure_id_path);
% mask and Kernel regression map for SSp
areaName = {'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n',...
    'SSp-bfd','SSp-un','VIS','RSP','VISa','VISrl','MO','ACA','PL'};
region_label = areaName(1:7);
% region_label = areaName(7);
st_region_indx = [];
for r = 1:numel(region_label)
    region_id = st.id(strcmp(st.acronym, region_label{r}));
    st_region_indx = [st_region_indx; find(cellfun(@(x)contains(x, sprintf('/%d/', region_id)), st.structure_id_path))];
end
[c,d] = ismember(projectedAtlas1,st_region_indx);
[row,col] = find(c);
ssp_index = [col,row];
%%
mn1{1} = 'AB_0004'; td1{1} = '2021-03-30'; en1{1} = 1;
mn1{2} = 'ZYE_0052'; td1{2} = '2021-12-18'; en1{2} = 2;
mn1{3} = 'ZYE_0056'; td1{3} = '2022-01-10'; en1{3} = 1;
mn1{4} = 'ZYE_0058'; td1{4} = '2022-03-09'; en1{4} = 1;
mn1{5} = 'ZYE_0012'; td1{5} = '2020-10-16'; en1{5} = 5;
%%
for kk = 1:5
    clear indx2 groupedCells filteredSpirals lia filteredSpirals_ssp U1 U2 meanTrace ...
        traceHilbert traceAmp traceAmp_mean spiral_amp spiral_amp_ssp
    mn = mn1{kk}; td = td1{kk}; en = en1{kk};
    formatOut = 'yyyymmdd';
    tdb = datestr(td,formatOut);
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    serverRoot = expPath(mn, td, en);
    wf_svd;
    %%
    folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetection4';
    dfolder = fullfile(folder,mn,td,num2str(en));
    load(fullfile(dfolder,'spiral_group.mat'));

    indx2 = cellfun(@(x) size(x,1)>3, archiveCell);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);

    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    [filteredSpirals(:,1),filteredSpirals(:,2)] = transformPointsForward(tform,filteredSpirals(:,1),filteredSpirals(:,2));
    filteredSpirals(:,1:2) = round(filteredSpirals(:,1:2));
    %%
    [lia,locb] = ismember(filteredSpirals(:,1:2),ssp_index,'rows');
    filteredSpirals_ssp = filteredSpirals(lia,:);
    %%
    scale = 8;
    U1 = U(1:scale:end,1:scale:end,1:50);
    mimg1 = mimg(1:scale:end,1:scale:end);
    U1 = U1./mimg1;
    U2 = reshape(U1,size(U1,1)*size(U1,2),size(U1,3));
    meanTrace = double(U2*dV(1:50,:));
    % filter 2-8Hz
    Fs = 35;
    meanTrace = meanTrace -mean(meanTrace ,2);
    % filter and hilbert transform work on each column
    meanTrace = meanTrace';
    [f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
    meanTrace = filtfilt(f1,f2,meanTrace);
    traceHilbert =hilbert(meanTrace);
    traceAmp = abs(traceHilbert);
    traceAmp = reshape(traceAmp,size(traceAmp,1),size(U1,1),size(U1,2));
    traceAmp_mean = mean(traceAmp,[2,3]);
    spiral_amp = traceAmp_mean(filteredSpirals(:,5));
    spiral_amp_ssp = traceAmp_mean(filteredSpirals_ssp(:,5));
    %%
    count = 1;
    for i = 10:10:100
        clear spiral_ampi
        spiral_ampi = spiral_amp(filteredSpirals(:,3)==i);
        mean_spiral_amp(count) = mean(spiral_ampi);
        std_spiral_amp(count) = std(spiral_ampi)/sqrt(numel(spiral_ampi));

        spiral_ampi_ssp = spiral_amp_ssp(filteredSpirals_ssp(:,3)==i);
        mean_spiral_amp_ssp(count) = mean(spiral_ampi_ssp);
        std_spiral_amp_ssp(count) = std(spiral_ampi_ssp)/sqrt(numel(spiral_ampi_ssp));

        count = count+1;
    end
    %%
    filteredSpirals(:,end+1) = spiral_amp;
    filteredSpirals_ssp(:,end+1) = spiral_amp_ssp;
    %%
    frameN = numel(t);
    h1 = figure('Renderer', 'painters', 'Position', [100 100 400 600]);
    subplot(2,1,1);
    shadedErrorBar(10:10:100,mean_spiral_amp*100,std_spiral_amp*100,'lineProps','b');
    hold on;
    shadedErrorBar(10:10:100,mean_spiral_amp_ssp*100,std_spiral_amp_ssp*100,'lineProps','g');
    xlabel('radius (pixels)');
    ylabel({'3-6 Hz amplitude','(df/f, %)'});

    subplot(2,1,2);
    color2 = cbrewer2('qual','Dark2',8);
    edges = [10:10:110];
    [N,edges] = histcounts(filteredSpirals(:,3),edges);
    N1 = N/frameN;
    bar(edges(1:end-1),N1,'FaceColor',color2(3,:))
    edges = [10:10:110];
    [N_ssp,edges] = histcounts(filteredSpirals_ssp(:,3),edges);
    N1_ssp = N_ssp/frameN;
    hold on;
    bar(edges(1:end-1),N1_ssp,'FaceColor',color2(5,:))
    xlabel('radius (pixels)');
    ylabel('spiral counts/ total frames');

    print(h1, [fname '_amp'], '-dpdf', '-bestfit', '-painters');
    save([fname '_amp'],'filteredSpirals','filteredSpirals_ssp','frameN','mean_spiral_amp','std_spiral_amp',...
        'mean_spiral_amp_ssp','std_spiral_amp_ssp','edges','N','N_ssp');
end
%% 
filteredSpirals_all = [];
filteredSpirals_ssp_all = [];
for kk = 1:5
    mn = mn1{kk}; td = td1{kk}; en = en1{kk};
    formatOut = 'yyyymmdd';
    tdb = datestr(td,formatOut);
    fname = [mn '_' tdb '_' num2str(en)];
    load([fname '_amp']);
    spiral_amp_all(kk,:) = mean_spiral_amp;
    spiral_amp_ssp_all(kk,:) = mean_spiral_amp_ssp;
    N_all(kk,:) = N/frameN;
    N_ssp_all(kk,:) = N_ssp/frameN;
    frameN_all(kk) = frameN;
    filteredSpirals(:,5) = filteredSpirals(:,5)+300000;
    filteredSpirals_ssp(:,5) = filteredSpirals_ssp(:,5)+300000;
    filteredSpirals_all = [filteredSpirals_all;filteredSpirals];
    filteredSpirals_ssp_all = [filteredSpirals_ssp_all;filteredSpirals_ssp];
end
%% seperated by session
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
subplot(1,2,1);
for i =1:5
    plot(10:10:100,spiral_amp_all(i,:)*100,'color',[0.5,0.5,0.5])
    hold on;
    plot(10:10:100,spiral_amp_ssp_all(i,:)*100,'color',[0 0 0])
    hold on;
end
xlabel('radius (pixels)');
ylabel({'3-6 Hz amplitude','(df/f, %)'});

subplot(1,2,2);
color2 = cbrewer2('qual','Dark2',8);
edges = [10:10:110];
N_all_mean = mean(N_all,1);
N_all_sem = std(N_all,1)/sqrt(5);

N_ssp_all_mean = mean(N_ssp_all,1);
N_ssp_all_sem = std(N_ssp_all,1)/sqrt(5);

bar(edges(1:end-1),N_all_mean,'FaceColor',color2(3,:))
hold on
er = errorbar(edges(1:end-1),N_all_mean,N_all_sem,N_all_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold on

bar(edges(1:end-1),N_ssp_all_mean,'FaceColor',color2(5,:))
hold on
er = errorbar(edges(1:end-1),N_ssp_all_mean,N_ssp_all_sem,N_ssp_all_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

xlabel('radius (pixels)');
ylabel('spiral counts/ total frames');

print(h1, 'spiral_radius_amplitude_by_session', '-dpdf', '-bestfit', '-painters');
%% combine all sessions and get amp histogram
count = 1;
for i = 10:10:100
    clear spiral_ampi
    spiral_ampi = filteredSpirals_all(filteredSpirals_all(:,3)==i,6);
    mean_spiral_amp(count) = mean(spiral_ampi);
    std_spiral_amp(count) = std(spiral_ampi)/sqrt(numel(spiral_ampi));

    spiral_ampi_ssp = filteredSpirals_ssp_all(filteredSpirals_ssp_all(:,3)==i,6);
    mean_spiral_amp_ssp(count) = mean(spiral_ampi_ssp);
    std_spiral_amp_ssp(count) = std(spiral_ampi_ssp)/sqrt(numel(spiral_ampi_ssp));

    count = count+1;
end
%%

h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
subplot(1,2,1);
shadedErrorBar(10:10:100,mean_spiral_amp*100,std_spiral_amp*100,'lineProps','b');
hold on;
shadedErrorBar(10:10:100,mean_spiral_amp_ssp*100,std_spiral_amp_ssp*100,'lineProps','g');
xlabel('radius (pixels)');
ylabel({'3-6 Hz amplitude','(df/f, %)'});

subplot(1,2,2);
color2 = cbrewer2('qual','Dark2',8);
edges = [10:10:110];
N_all_mean = mean(N_all,1);
N_all_sem = std(N_all,1)/sqrt(5);

N_ssp_all_mean = mean(N_ssp_all,1);
N_ssp_all_sem = std(N_ssp_all,1)/sqrt(5);

bar(edges(1:end-1),N_all_mean,'FaceColor',color2(3,:))
hold on
er = errorbar(edges(1:end-1),N_all_mean,N_all_sem,N_all_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold on

bar(edges(1:end-1),N_ssp_all_mean,'FaceColor',color2(5,:))
hold on
er = errorbar(edges(1:end-1),N_ssp_all_mean,N_ssp_all_sem,N_ssp_all_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

xlabel('radius (pixels)');
ylabel('spiral counts/ total frames');
print(h1, 'spiral_radius_amplitude_all', '-dpdf', '-bestfit', '-painters');