function spiral_nomatch = getNoMatchedSprials(spirals_raw,spirals_predict)
% only care about the situation when only one spiral is within this
% distance
spiral_nomatch = [];                                                % check left hemisphere first 
for i = 1:size(spirals_predict)
    clear spiral_p1 index1 spiral_index
    frame = spirals_predict(i,5);
    spiralx = spirals_predict(i,1); 
    spiraly = spirals_predict(i,2);

    spiral_index = find(spirals_raw(:,5) == frame);
    spiral_p1 = spirals_raw(spiral_index,:);
    index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
        spiral_p1(:,1)<= spiralx+100 & ...
        spiral_p1(:,2)>= spiraly-100 & ...
        spiral_p1(:,2)<= spiraly+100);
    if isempty(index1)
        spiral_temp = [spirals_predict(i,1:5)];
        spiral_nomatch = [spiral_nomatch;spiral_temp];
    end
end