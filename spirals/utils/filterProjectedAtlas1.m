function [projectedAtlas1,projectedTemplate1] = filterProjectedAtlas1(projectedAtlas,projectedTemplate,st)
spath = string(st.structure_id_path);
% only select isocortex in the atlas
projectedAtlas1 = projectedAtlas;
projectedTemplate1 = projectedTemplate;
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/688/695/315/'); 
% spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;