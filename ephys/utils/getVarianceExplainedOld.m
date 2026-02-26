%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new.csv'));
%%
T1 = T(30:35,:);
%%
scale = 4;
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\variance_summary\MB';
for kk = 1:6
    %%
    clearvars -except T1 kk data_folder save_folder folder scale
    %% load widefield and spiking data
    ops = get_session_info2(T1,kk,data_folder);
    [U,V,t,mimg] = loadUVt1(ops.session_root);                             % load U,V, t
    Ut = U(1:scale:end,1:scale:end,1:50);
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    [sp] = loadKSdir2(ops.session_root);
    [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t);
    WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
    dV1 = double(dV(:,~isnan(WF2ephysT1)));
    V1 = double(V(1:50,~isnan(WF2ephysT1)));
    t1 = t(~isnan(WF2ephysT1));
    [MUA_std] = get_MUA_bin(sp,WF2ephysT);
    dV1 = double(dV1(1:50,:));
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    %%
    filename = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_explained_var'];
    ffname = fullfile(folder,filename);
    %% full prediction
    [dV_raw,dV_predict,epoch_indx] = get_prediction(dV1,MUA_std);
    explained_var_all = get_variance_explained(Ut,dV_raw,dV_predict);
    save(ffname,'explained_var_all');
end