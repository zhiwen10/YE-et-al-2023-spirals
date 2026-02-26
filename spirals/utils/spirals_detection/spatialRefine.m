function [spiralCT1] = spatialRefine(tracePhase,pwframe,params)
pgridx = params.pgridx;
pgridy = params.pgridy;
spiralRange = params.spiralRange;
th = params.th;
spiralCT1 = [];
kcount = 1;

for ixy = 1:size(pwframe,1)
    px = pwframe(ixy,1); py = pwframe(ixy,2);
    roi1 = tracePhase((py-20):(py+20),(px-20):(px+20));
    %% search for spirals within 20x20 pixels fine grid
    % only need to conform by a single radius
    % then average across all candidate centers for true spiral center
    rs1 = 2;
    spiralD = zeros(numel(pgridx(:)),1);
    for kkk = 1:numel(pgridx(:))
        px1 = pgridx(kkk); py1 = pgridy(kkk);
        ii = 1;
        while spiralD(kkk) == 0 & ii<2
            r = rs1(ii);
            cx = round(r*cosd(th)+px1);
            cy = round(r*sind(th)+py1); 
            %%
            ph = roi1(sub2ind(size(roi1),cy, cx));
            phdiff = angdiff(ph);
            ph2(1) = ph(1);
            for i = 2:numel(ph)                
                ph2(i) = [ph2(i-1)+phdiff(i-1)];
            end 
            ph3 = abs(ph2-ph2(1));                
            [N,edges] = histcounts(ph,spiralRange);
            AngleRange = ph3(end);       
            if AngleRange>5 & AngleRange<7 & all(N)
                spiralD(kkk) = 1;
            end
            ii = ii+1;
        end
    end
    %%
    pgridx1 = pgridx(logical(spiralD));
    pgridy1 = pgridy(logical(spiralD));
    spiralCT = [];
    if sum(spiralD)>0
        spiralCT = round(mean([pgridx1,pgridy1],1));
        spiralCT1(kcount,1) = px-21+spiralCT(1);
        spiralCT1(kcount,2) = py-21+spiralCT(2);
        kcount = kcount+1;
    end
end