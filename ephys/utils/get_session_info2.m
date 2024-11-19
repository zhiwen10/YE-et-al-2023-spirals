function ops = get_session_info2(T,kk,data_folder)
%% session info
mn1 = T.MouseID;
td1 = T.date;
en1 = T.folder;
imec1 = T.imec;
probeName1 = T.probeName;
doubleLength1 = T.doubleLength;
%%
ops.mn = mn1{kk}; ops.tda = td1(kk); ops.en = en1(kk); 
ops.imec = imec1(kk); ops.probeName = probeName1{kk};
ops.doubleLength = doubleLength1(kk);
imecN = 1;
ops.td = datestr(ops.tda,'yyyy-mm-dd');
ops.tdb = datestr(ops.tda,'yyyymmdd');
ops.fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
%%
subfolder = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
ops.session_root = fullfile(data_folder,'ephys','svd_spikes',subfolder);
%% load data parameters
if ops.doubleLength == 1
    ops.chanMap = fullfile(data_folder,...
        'ephys','config_files','NPtype24_doubleLengthStripe_botRow0_ref0.mat');     % NP2.4 doubleLength
elseif ops.doubleLength == 0
    ops.chanMap = fullfile(data_folder,...
        'ephys','config_files','NPtype24_quadrupleLengthStripe_botRow0_ref0.mat');  % NP2.4 quadrupleLength
elseif ops.doubleLength == 2
    ops.chanMap = fullfile(data_folder,...
        'ephys','config_files','neuropixPhase3B1_kilosortChanMap.mat');
elseif ops.doubleLength == 3
    ops.chanMap = fullfile(data_folder,...
        'ephys','config_files','NPtype24_hStripe_botRow0_ref1.mat');
elseif ops.doubleLength == 4
    ops.chanMap = fullfile(data_folder,...
        'ephys','config_files','NPtype21_botRow0.mat');
end
%%
ops.chanMap1 = chanMapReorder(ops.chanMap);
ops.curation = T.Curation(kk);
end
