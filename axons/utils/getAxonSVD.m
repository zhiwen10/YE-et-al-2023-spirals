function [soma_center,axon_vector, angle1, polarity] = getAxonSVD(axon_terminal_all,soma_all)
%% function to get firt PC of axon branching points
% [axon_terminal_all] is cell array of 3d axon coordinates;
% [soma_all] is cell array of 3d soma coordinates;

% [vector_all] is the eigenvector of first PC  
% [angle1] is the angle of first PC
% [polarity] is the ratio of eigenvalue of first and second PC -1
% so that the minumum of polarity is zero, because by default first PC
% value is >= second PC.

%%
scale = 1;
% a = cellfun(@(x) not(isempty(x)), axon_terminal_all);
% axon_terminal_all = axon_terminal_all(a);
% soma_all = soma_all(a);
cell_n = numel(axon_terminal_all);  
for icell = 1:cell_n
    axon_current = axon_terminal_all{icell};
    axon_current_2d = double(axon_current(:,[3,1]));
    
    soma_current_2d = soma_all{icell};
    soma_current_2d = soma_current_2d(1,:);
    soma_center(icell,:) = double(soma_current_2d(:,[3,1]));
    
    % zero mean beofre svd
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
    
    % axon svd
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag(:,icell) = diag(S);
    axon_vector(icell,:) = [U(1,1),U(2,1)];
    angle1(icell,1) = round(rad2deg(atan(U(2,1)/U(1,1))))+1;
    angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
    % plot(soma_center(1),soma_center(2),'.','MarkerFaceColor',color2(angle1,:),'MarkerSize',16)
end
polarity = S_diag(1,:)./S_diag(2,:)-1;
% polarity = S_diag(1,:);
if ~iscolumn(polarity)
    polarity = polarity';
end
% polarity = ones(size(vector_all,1),1);