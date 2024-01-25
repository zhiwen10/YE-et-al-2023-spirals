githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
% addpath(genpath('C:\Users\Steinmetz lab\Documents\git\AP_scripts_cortexlab-master'))
% addpath(genpath(fullfile(githubDir, 'wheelAnalysis')))
% addpath(genpath(fullfile(githubDir, 'PatternDetection')))
% addpath(genpath(fullfile(githubDir, 'spiralDetection')));
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\shadedErrorBar'));
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
T = readtable('session_list_sorted2.csv');
T1 = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
%%
prediction_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\prediction';
save_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow_new1\batch_flow_code';
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\roi';
%%
for kk = 1:size(T1,1) 
    %%
    ops = get_session_info(T1,kk);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    fname_predict = fullfile(prediction_folder,[fname '_dv_predict.mat']);
    load(fname_predict,'explained_var_all');
    load(fullfile(roi_folder,[fname '_roi.mat']));
    [Ut,mimg,V1,dV1,MUA_std] =get_wf_mua1(ops);     
    scale = 4;
    mimg1 = mimg(1:scale:end,1:scale:end);
    %%
    ffname = fullfile(save_folder,[fname '.mat']); % save folder    
    explained_var_all(explained_var_all<0) = 0;
    var_mean = squeeze(mean(explained_var_all,3));
    BW1 = (var_mean>=0.1);
    bw = poly2mask(roi.Position(:,1),roi.Position(:,2),size(mimg,1)+240,size(mimg,2)+240);
    bw = bw(121:size(bw,1)-120,121:size(bw,2)-120);
    bw = bw(1:scale:end,1:scale:end);
    BW = (BW1&bw);
    %%
    len = 5000; % sample size to use
    [N,edges,phase_mu,phase_var,flow_mu,flow_var,traceAmp_mean] = compare_flow1(Ut,mimg1,dV1,MUA_std,BW,len);   
    save(ffname,'N','edges','phase_mu','phase_var','flow_mu','flow_var','traceAmp_mean');
end