%%
phasemap1 = readNPY('C:\Users\Steinmetz lab\OneDrive - UW\Documents\git\Ye-et-al-2023-spirals\simulation\mirrored\mirrored_symmetry.npy');
phasemap2 = readNPY('C:\Users\Steinmetz lab\OneDrive - UW\Documents\git\Ye-et-al-2023-spirals\simulation\mirrored\mirrored_symmetry_control.npy');
%%
phasemap1a = phasemap1(:,:,1:2:20);
phasemap2a = phasemap2(:,:,1:2:20);
h1b = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);
for i = 1:10
    ax1 = subplottight(2,10,i);
    ax1.Position(1) = ax1.Position(1)+0.005;
    ax1.Position(2) = ax1.Position(2);
    ax1.Position(3) = ax1.Position(3)-0.01;
    ax1.Position(4) = ax1.Position(4)-0.01;
    phasetemp = squeeze(phasemap2a(:,:,i));
    phasetemp = flipud(phasetemp);
    imagesc(phasetemp);
    colormap(ax1,colorcet('C06'));
    hold on;
    plot([100,100],[0,200],'--w');
    hold on;
    plot([0,200],[67,67],'--w');
    axis image; axis off;
end

for i = 1:10
    ax1 = subplottight(2,10,10+i);
    ax1.Position(1) = ax1.Position(1)+0.005;
    ax1.Position(2) = ax1.Position(2);
    ax1.Position(3) = ax1.Position(3)-0.01;
    ax1.Position(4) = ax1.Position(4)-0.01;
    phasetemp = squeeze(phasemap1a(:,:,i));
    phasetemp = flipud(phasetemp);
    imagesc(phasetemp);
    colormap(ax1,colorcet('C06'));
    hold on;
    plot([100,100],[0,200],'--w');
    hold on;
    plot([0,200],[67,67],'--w');
    axis image; axis off;
end
%%
print(h1b, 'mirrored_spirals_model.pdf','-dpdf', '-bestfit', '-painters');