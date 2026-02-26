function [pwAll5] = spiralRadiusCheck2(tracePhase,pwframe,params)
th = params.th;
spiralRange = params.spiralRange;
rsRCheck = params.rsRCheck;
%%
pwAll5 = [];
for ixy = 1:size(pwframe,1)
    px = pwframe(ixy,1); py = pwframe(ixy,2);
    % rs1 = 10:10:100;
    spiralR = zeros(numel(rsRCheck),2);
                %%
            for ii = 1:numel(rsRCheck)
                    r = rsRCheck(ii);
                    cx = round(r*cosd(th)+px);
                    cy = round(r*sind(th)+py); 
                    %%
                    ph = tracePhase(sub2ind(size(tracePhase),cy, cx));  
                    if not(any(isnan(ph)))
                        phdiff = angdiff(ph);
                    else
                        indx1 = isnan(ph);
                        ph(indx1) = phAll(ii-1,indx1);
                    end
                    phAll(ii,:) = ph;
                    ph2(1) = ph(1);
                    for i = 2:numel(ph)                
                        ph2(i) = [ph2(i-1)+phdiff(i-1)];
                    end 
                    ph2 = ph2-ph2(1);
                    ph2All(ii,:) = ph2(end);
                    ph3 = abs(ph2);                
                    [N,edges] = histcounts(ph,spiralRange);
                    AngleRange = ph3(end);       
                    if AngleRange>4 & AngleRange<8 & all(N)
                        spiralR(ii,1) = 1;
                        if ph2(end)>0
                            spiralR(ii,2) = 1; % counterclock wise
                        else
                            spiralR(ii,2) = -1;
                        end
                    end
            end
            if sum(spiralR(:,1),1)>0
                % maxR = rsRCheck(find(spiralR(:,1),1,'last'));
                indxZ = find(spiralR(:,1),1,'first');
                spiralR1 = spiralR;
                if indxZ>1                  
                    spiralR1(1:indxZ-1,2) = spiralR1(indxZ,2);
                end
%                 indxZ1 = find(not(spiralR1(:,1)),1,'first');
%                 maxR = rsRCheck(indxZ1-1);
%                 direction = spiralR(indxZ1-1,2);

                % find raduis where it changed direction as end of spiral
                % otherwise is max radius of 100
                indxZ1 = find(spiralR1(:,2)~=spiralR1(indxZ,2),1,'first'); 
                if not(isempty(indxZ1))
                    maxR = rsRCheck(indxZ1-1);
                else
                    maxR = rsRCheck(end);
                end
                % direction is set as the direction at the smallest radius 
                direction = spiralR(indxZ,2);
                pwTemp = [px,py,maxR,direction];
                pwAll5 = [pwAll5;pwTemp];
            end
        end
end