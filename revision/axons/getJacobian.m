function [itype] = getJacobian(vx,vy)
dxx = mean(vx(:,2) - vx(:,1));
dxy = mean(vx(1,:) - vx(2,:));
dyx = mean(vy(:,2) - vy(:,1));
dyy = mean(vy(1,:) - vy(2,:));
jacobians = [dxx dxy; dyx dyy];
det1 = det(jacobians);
trace1 = trace(jacobians);
% Classify critical point by its Jacobian
if det(jacobians) < 0
    itype = 'saddle';
elseif det(jacobians) > 0 & trace(jacobians)^2 >= 4*det(jacobians)
    itype = 'nodes';
elseif det(jacobians) > 0 & trace(jacobians)^2 < 4*det(jacobians)
    itype = 'spiral';
% elseif det(ijac) > 0 & trace(ijac)^2 == 4*det(ijac)
%     itype = 'stars';
end