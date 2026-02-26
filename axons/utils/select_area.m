function [index,Uselected] =  select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale)
%% input area atlas path, find and select them in atlas 
% areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
% areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
% areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
% areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
% areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
% areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
% areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
% areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
% areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
% areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
% sensoryArea = strcat(areaPath(:));

% st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
% spath = string(st.structure_id_path);

% projectedAtlas2 = projectedAtlas1;
% projectedTemplate2 = projectedTemplate1;
 
%%
if not(isempty(sensoryArea))
    spath3 = startsWith(spath,sensoryArea);
    idFilt1 = st.index(spath3);
    [Lia1,Locb1] = ismember(projectedAtlas1,idFilt1);
    projectedAtlas1(~Lia1) = 0; 
    projectedTemplate1(~Lia1) = 0;
end

if strcmp(hemi,'right')
    projectedAtlas1(:,1:size(projectedAtlas1,2)/2) = 0; 
    projectedTemplate1(:,1:size(projectedTemplate1,2)/2) = 0;
elseif strcmp(hemi,'left')
    projectedAtlas1(:,size(projectedAtlas1,2)/2+1:end) = 0; 
    projectedTemplate1(:,size(projectedAtlas1,2)/2+1:end) = 0; 
else
    disp('');  % do nothing
end
scale1 = 1;
% plotOverlaidAtlas(projectedAtlas1,projectedTemplate1,coords,scale1);
projectedAtlas2 = projectedAtlas1(1:scale:end,1:scale:end);
index = find(projectedAtlas2);
UregDown2d = Utransformed(1:scale:end,1:scale:end,:);
UregDown1d = reshape(UregDown2d,size(UregDown2d,1)*size(UregDown2d,2),size(UregDown2d,3));
Uselected = UregDown1d(index,:);