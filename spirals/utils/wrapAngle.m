function ph2= wrapAngle(ph)
phdiff = angdiff(ph);
ph2(1) = ph(1);
for i = 2:numel(ph)                
    ph2(i) = [ph2(i-1)+phdiff(i-1)];
end 
end