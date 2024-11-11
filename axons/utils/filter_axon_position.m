function [axonCoords_filt] = filter_axon_position(axonCoords,av,st_region_indx,hemi)
axon_area = [];
axonCoords2 = round(axonCoords(:,2:4)/10); % compatibility with av (10 micron resolution)
% fixing rounding errors:
axonCoords2(axonCoords2<1) = 1;
axonCoords2(axonCoords2(:,1)>size(av,1),1) = size(av,1); axonCoords2(axonCoords2(:,2)>size(av,2),2) = size(av,2); 
axonCoords2(axonCoords2(:,3)>size(av,3),3) = size(av,3); 
for nIdx = 1:size(axonCoords,1)
    currentAxon = axonCoords2(nIdx,:);
    avIdxAxon = sub2ind(size(av), currentAxon(1), currentAxon(2), currentAxon(3));
    stIdxAxon = av(avIdxAxon);
    axon_area = [axon_area;stIdxAxon];
end
%% only use left hemisphere
if hemi== [] % both hemisphere
    axonCoords_filt = axonCoords2(ismember(axon_area, st_region_indx),:);
elseif hemi==1 % left hemisphere only
    axonCoords_filt = axonCoords2(ismember(axon_area, st_region_indx) & axonCoords2(:,3)<size(av,3)/2,:);
elseif hemi == 2 % right hemisphere only 
    axonCoords_filt = axonCoords2(ismember(axon_area, st_region_indx) & axonCoords2(:,3)>size(av,3)/2,:);
end
