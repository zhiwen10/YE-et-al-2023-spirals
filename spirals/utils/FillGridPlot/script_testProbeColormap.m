
%%

chanMap = load('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref0.mat');


%%

% your data should be size [nFrames nChannels] where the channels are in
% the same order that the coordinates are in the channel map. so don't
% resize the data to be 8 x 48
data = rand(20, numel(chanMap.xcoords)); 
cax = [0 1];
colMap = hsv(100);
%%
f = figure(1); 
ax = gca;  
f.Color = 'k';

ps = makeProbePlot(ax, chanMap, 12);
axis image;
% axis off; 
%%
for d = 1:size(data,1)
    colorSites(ps, data(d,:), colMap, cax);
    drawnow;
end
