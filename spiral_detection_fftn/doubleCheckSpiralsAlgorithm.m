function [pwAll] = doubleCheckSpiralsAlgorithm(tracePhase,pwframe,params)
% A = squeeze(tracePhase(frame,:,:));
% frameN = frame+frameStart-1;
%%
rs = params.rs;
th = params.th;
spiralRange = params.spiralRange;
pwAll = [];
%%
for ixy = 1:size(pwframe,1)
    px = pwframe(ixy,1); py = pwframe(ixy,2);
    pw = [];
    pwTemp = [];
    ph2 = zeros(1,10);
    ph3 = zeros(1,10);
    for rn = 1:numel(rs)
        r = rs(rn);
        [spiralTemp] = checkSpiral(tracePhase,px,py,r,th,spiralRange);
        pw = [pw;spiralTemp];
    end   
    % if 2 out of 3 radius satisfy the criteria, then a candidate spiral center
    if sum(pw(:,4))>=2
            pwTemp1 = [px,py];   
            pwAll = [pwAll;pwTemp1];
    end   
end