function [sp] = get_spikes(ksSubFolder,manulacluster)
sp = loadKSdir(ksSubFolder);
%% good cluster only cgs  == 2
if nargin<2
    manulacluster = 1;
end
if manulacluster ==1
    % if manually curated in phy2
    sp.gcluster  = sp.cids(sp.cgs==2);
else
    load([ksSubFolder '\rez.mat']);
    sp.gcluster = find(rez.good==1)-1;
end

[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
% driftmap
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksSubFolder);
sp.spikeTimes = spikeTimes;
sp.spikeAmps = spikeAmps;
sp.spikeDepths = spikeDepths;
sp.spikeSites = spikeSites;
end