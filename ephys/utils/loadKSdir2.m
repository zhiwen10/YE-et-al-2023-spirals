function spikeStruct = loadKSdir2(ksDir, varargin)
if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end
if ~isfield(params, 'excludeNoise')
    params.excludeNoise = true;
end
if ~isfield(params, 'loadPCs')
    params.loadPCs = false;
end
% load spike data
ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/30000;                                                     % sampling rate: 30k Hz
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy'));          % note: zero-indexed
clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
temps = readNPY(fullfile(ksDir, 'templates.npy'));
winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

spikeTimes = readNPY(fullfile(ksDir, 'spikeTimes.npy'));
spikeAmps = readNPY(fullfile(ksDir, 'spikeAmps.npy'));
spikeDepths = readNPY(fullfile(ksDir, 'spikeDepths.npy'));
spikeSites = readNPY(fullfile(ksDir, 'spikeSites.npy'));

cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv')) 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
elseif exist(fullfile(ksDir, 'cluster_group.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
elseif exist(fullfile(ksDir, 'cluster_KSLabel.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_KSLabel.tsv');
end 

if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);
    
    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));     
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));       
    end
end

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);

spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.spikeAmps = spikeAmps;
spikeStruct.spikeDepths = spikeDepths;
spikeStruct.spikeSites = spikeSites;
spikeStruct.spikeTimes = spikeTimes;
spikeStruct.gcluster  = cids(cgs==2);