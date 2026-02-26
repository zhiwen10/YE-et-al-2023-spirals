function getPhaseFlowMatchingIndex(T,data_folder,save_folder)
%%
for kk = 1:size(T,1) 
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    fname_predict = fullfile(data_folder,'ephys','dv_prediction',...
        [fname '_dv_predict.mat']);
    load(fname_predict,'explained_var_all');
    load(fullfile(data_folder,'ephys','roi',[fname '_roi.mat']));
    [Ut,mimg,V1,dV1,MUA_std] =get_wf_mua2(ops);     
    scale = 4;
    mimg1 = mimg(1:scale:end,1:scale:end);
    %% 
    explained_var_all(explained_var_all<0) = 0;
    var_mean = squeeze(mean(explained_var_all,3));
    BW1 = (var_mean>=0.1);
    bw = poly2mask(roi.Position(:,1),roi.Position(:,2),...
        size(mimg,1)+240,size(mimg,2)+240);
    bw = bw(121:size(bw,1)-120,121:size(bw,2)-120);
    bw = bw(1:scale:end,1:scale:end);
    BW = (BW1&bw);
    %%
    len = 5000;                                                            % sample size to use
    [N,edges,phase_mu,phase_var,flow_mu,flow_var,traceAmp_mean] = ...
        compare_flow1(Ut,mimg1,dV1,MUA_std,BW,len);   
    %% save data
    ffname = fullfile(save_folder,[fname '.mat']);                         %
    save folder   
    save(ffname,'N','edges','phase_mu','phase_var',...
        'flow_mu','flow_var','traceAmp_mean');
end