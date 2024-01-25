function [MUA_std] = get_MUA_bin(sp,WF2ephysT)
%% bin MUA based on WF frame time
for m = 1:numel(sp.gcluster)
    clear incl 
    % incl = (sp.spikeAmps>20 &(sp.clu==sp.gcluster(m)) & sp.spikeDepths<depth);
    incl = (sp.spikeAmps>20&(sp.clu==sp.gcluster(m)));
    spikeTimes2 = sp.spikeTimes(incl); 
    spikeT{m} = spikeTimes2;
end
%% 
for m = 1:numel(sp.gcluster)
    SpikeSite = spikeT{m};
    [n] = histcounts(SpikeSite, WF2ephysT);
    n = [0 n];
    MUA(:,m) = n; 
end
%%
 MUA = reshape(MUA,size(MUA,1),size(MUA,2)*size(MUA,3));
 MUA = MUA';
 MUA_std = (MUA-mean(MUA,2))./nanstd(MUA,[],2);
 MUA_std = MUA_std(std(MUA_std,[],2)>0,:);
 % MUA_std(end+1,:) = 1:size(MUA_std,2); %add time shift variant
end