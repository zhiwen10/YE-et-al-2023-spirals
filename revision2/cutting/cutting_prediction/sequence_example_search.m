%%
spiral_length3 = cellfun(@(x) size(x,1), archiveCell);
archiveCell2 = archiveCell(spiral_length3==10);
seq1 = archiveCell2{2};
frames = seq1(1,5):seq1(end,5);
centers2 = centers/8; radius2 = radius/8;
[h1b,r_var] = plotSpiralSequences2(Upt,dV_predict,U1t,dV1,t1,centers2,radius2,frames);
%%
figure;
BW2 = BW(1:8:end,1:8:end);
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/8;
th2 = 1:5:360; 
im_phase = imagesc(1-r_var);
axis image; axis off;
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([0,1]);

for kk = 1:size(centers,1)
    px1 = centers2(kk,1);
    py1 = centers2(kk,2);
    r = radius2(kk,1);
    cx2 = round(r*cosd(th2)+px1);
    cy2 = round(r*sind(th2)+py1);
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end