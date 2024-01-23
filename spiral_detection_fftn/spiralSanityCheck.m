%%
function spiralSanityCheck(U,dV,t,pwAll,frame2,framei,params,rate1)
th = [1:36:360]; 
%%
dV1 = dV(1:50,frame2);
t1 = t(frame2);
[trace2d2,traceAmp2,tracePhase2] = spiralPhaseMap4(U(:,:,1:50),dV1,t1,params,rate1);
figure;
imagesc(squeeze(tracePhase2(:,:,framei)));
colormap(hsv);
axis image;
hold on; 
framei = frame2(framei);
pwi = pwAll(pwAll(:,end)==framei,:);
pwi = pwi*params.downscale;
hold on;
scatter(pwi(:,1),pwi(:,2),32,'k','filled');
for i = 1:size(pwi,1)
    px = pwi(i,1);
    py = pwi(i,2);
    cx = round(pwi(i,3)*cosd(th)+px);
    cy = round(pwi(i,3)*sind(th)+py); 
    hold on;
    if pwi(i,4)==1
        c = 'k';
    else
        c = 'r';
    end
    scatter(cx,cy,16,c,'filled')
end
