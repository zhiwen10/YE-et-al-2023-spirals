function spiral_count_sum = sortSpirals(pwAll,indexSSp,frames_stimOn)
[lia2,locb2] = ismember(pwAll(:,1:2),indexSSp,'rows');
pwAll2 = pwAll(lia2,:);
spiral_cell2 = cell(size(frames_stimOn));
for j = 1:size(frames_stimOn,1)
    for k = 1:size(frames_stimOn,2)
        clear indx1 current_spirals
        current_frame = frames_stimOn(j,k);
        indx1 = find(pwAll2(:,5) == current_frame);
        current_spirals = pwAll2(indx1,:);
        spiral_cell2{j,k} = current_spirals;
    end
end  
spiral_count_sum = SpiralCountBins(spiral_cell2); 