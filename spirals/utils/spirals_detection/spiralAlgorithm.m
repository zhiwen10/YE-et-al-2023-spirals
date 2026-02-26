function [pwAll1] = spiralAlgorithm(tracePhase,params)
xxRoi =params.xxRoi;
yyRoi = params.yyRoi;
rs = params.rs;
th = params.th;
spiralRange = params.spiralRange;
dThreshold = params.dThreshold;
%%
clear ph phdiff 
pwAll1 = [];
for ixy = 1:numel(xxRoi)
        px = xxRoi(ixy); py = yyRoi(ixy);
        pw = [];
        pwTemp = [];
        ph2 = zeros(1,10);
        ph3 = zeros(1,10);
        for rn = 1:numel(rs)
            r = rs(rn);
            [spiralTemp] = checkSpiral(tracePhase,px,py,r,th,spiralRange);
            pw = [pw;spiralTemp];
        end
        %% if 2 out of 3 radius satisfy the criteria, then a candidate spiral center                       
        if sum(pw(:,4))>=2
                % frameN = frame+frameStart-1;
                pwTemp1 = [px,py];   
                pwAll1 = [pwAll1;pwTemp1];
        end                
end
end