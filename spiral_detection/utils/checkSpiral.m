function [spiralTemp] = checkSpiral(A,px,py,r,th,spiralRange)
% th = 1:36:360; 
% spiralRange = linspace(-pi,pi,5);
cx = round(r*cosd(th)+px);
cy = round(r*sind(th)+py);
%%
ph = A(sub2ind(size(A),cy, cx));
phdiff = angdiff(ph);
ph2(1) = ph(1);
for i = 2:numel(ph)                
    ph2(i) = [ph2(i-1)+phdiff(i-1)];
end 
ph3 = abs(ph2-ph2(1));                
[N,edges] = histcounts(ph,spiralRange);
%%
AngleRange = abs(ph2(end)-ph2(1));
spiralTemp = [px,py,r,0];
if AngleRange>5 & AngleRange<7 & all(N)
    spiralTemp(1,4) = 1;                  
end
end