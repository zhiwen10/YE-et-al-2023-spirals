function ops = get_session_info(T,kk)
%% session info
mn1 = T.mouseId;
td1 = T.date;
en1 = T.en;
imec1 = T.imec;
probeName1 = T.probeName;
doubleLength1 = T.doubleLength;
%%
ops.mn = mn1{kk}; ops.tda = td1(kk); ops.en = en1(kk); ops.imec = imec1(kk); ops.probeName = probeName1{kk};
ops.doubleLength = doubleLength1(kk);
imecN = 1;
ops.td = datestr(ops.tda,'yyyy-mm-dd');
ops.tdb = datestr(ops.tda,'yyyymmdd');
ops.fname = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_inverse'];
%%
ops.serverRoot = expPath(ops.mn, ops.td, ops.en);
%%
ops.ksSubFolder = fullfile(fileparts(getProbeFile(ops.serverRoot, ops.probeName,ops.imec,1)));
%% load data parameters
if ops.doubleLength == 1
    ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_doubleLengthStripe_botRow0_ref0.mat'; % NP2.4 doubleLength
elseif ops.doubleLength == 0
    ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_quadrupleLengthStripe_botRow0_ref0.mat'; % NP2.4 quadrupleLength
elseif ops.doubleLength == 2
    ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\neuropixPhase3B1_kilosortChanMap.mat';
elseif ops.doubleLength == 3
    ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
elseif ops.doubleLength == 4
    ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype21_botRow0.mat';
end
%%
ops.chanMap1 = chanMapReorder(ops.chanMap);
%%
ops.curation = T.Curation(kk);
end
