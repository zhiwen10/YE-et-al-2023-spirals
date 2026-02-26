%% test gradient
% spiral
vx = [1,0;0,-1];
vy = [0,-1;1,0];
% saddle
% vx = [1,-1;1,-1];
% vy = [1,1,;-1,-1];
% sink/source star
% vx = [1,-1;1,-1];

% vy = [-1,-1;1,1];
% sink/source star
% vx = [-1,1;-1,1];
% vy = [1,1;-1,-1];
% sink/source node
% vx = [1,-1;1,-1]
% vy = [-2,-1;1,2];

dxx = mean(vx(:,2) - vx(:,1));
dxy = mean(vx(1,:) - vx(2,:));
dyx = mean(vy(:,2) - vy(:,1));
dyy = mean(vy(1,:) - vy(2,:));
jacobians = [dxx dxy; dyx dyy];
% Classify critical point by its Jacobian
if det(jacobians) < 0
    itype = 'saddle';
elseif det(jacobians) > 0 & trace(jacobians)^2 >= 4*det(jacobians)
    itype = 'sink/source';
elseif det(jacobians) > 0 & trace(jacobians)^2 < 4*det(jacobians)
    itype = 'spiral';
% elseif det(ijac) > 0 & trace(ijac)^2 == 4*det(ijac)
%     itype = 'stars';
end