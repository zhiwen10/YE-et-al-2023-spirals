function ax = plotAngleStairsShuffle(ax,edges,N1,N_rp)

mean_N_rp = mean(N_rp,1);
% std_N_rp = std(N_rp,[],1)/sqrt(size(N_rp,1));
std_N_rp = std(N_rp,[],1);
stairs(ax,edges(1:end-1),N1,'k');
% hold on;
% stairs(ax,edges(1:end-1),mean_N_rp+std_N_rp,'k');
% hold on;
% stairs(ax,edges(1:end-1),mean_N_rp-std_N_rp,'k');
hold on;
stairs(ax,edges(1:end-1),mean_N_rp,'color',[0.5,0.5,0.5]);

%Calculate the error bars
uE = mean_N_rp+std_N_rp;
lE = mean_N_rp-std_N_rp;
uEs = uE(1:end-1);
lEs = lE(1:end-1);
uE_new(1:2:2*numel(uE)-1) = uE;
uE_new(2:2:2*numel(uE)-1) = uEs;
lE_new(1:2:2*numel(lE)-1) = lE;
lE_new(2:2:2*numel(lE)-1) = lEs;

xu = edges(1:end-1);
xus = xu(2:end);
xu_new(1:2:2*numel(xu)-1) = xu;
xu_new(2:2:2*numel(xu)-1) = xus;

%Make the patch (the shaded error bar)
yP=[lE_new,fliplr(uE_new)];
xP=[xu_new,fliplr(xu_new)];
Hpatch=patch(ax,xP,yP,1);
patchColor=[0.5,0.5,0.5]; faceAlpha = 0.5;
set(Hpatch,'facecolor',patchColor, ...
    'edgecolor','none', ...
    'facealpha',faceAlpha, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')

% xticks(ax,[0 30 60 90 120 150 180])
% xticklabels(ax,{'0' '30' '60' '90' '120' '150' '180'})
xticks(ax,[0,45,90,135,180]);
xticklabels(ax,{'0','1/4*pi','1/2*pi','3/4*pi','pi'});

xlabel('angle bins');
ylabel('bin counts');
legend({'data',['shuffle' char(177) 'std']});