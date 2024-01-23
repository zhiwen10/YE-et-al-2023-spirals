%% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath(fullfile(githubdir2, 'allenCCF')));
addpath(genpath(fullfile(githubdir2, 'AP_scripts_cortexlab-master')));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
session_all = find(T.use);
session_total = numel(session_all);
%%
for kk = 1:15
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    load(fullfile(serverRoot, 'blue', ['dataSummary.mat']));
    cumsumSv = cumsum(Sv);
    Sv_ratio(:,kk) = cumsumSv/totalVar;
end
%%
mean_Sv_ratio = mean(Sv_ratio,2);
std_Sv_ratio = std(Sv_ratio,0,2);
sem_Sv_ratio = std_Sv_ratio/sqrt(size(Sv_ratio, 2));
h1= figure;
shadedErrorBar(1:500, mean_Sv_ratio(1:500), sem_Sv_ratio(1:500), 'lineprops', '-r')
xline(50,'--');
text(100,0.96,num2str(round(mean_Sv_ratio(50)*100)/100));
%%
print(h1, 'data_variance', '-dpdf', '-bestfit', '-painters');