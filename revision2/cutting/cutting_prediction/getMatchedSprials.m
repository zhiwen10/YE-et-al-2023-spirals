function spiral_match = getMatchedSprials(spirals_raw,spirals_predict)
% only care about the situation when only one spiral is within this
% distance
spiral_match = [];                                                % check left hemisphere first 
for i = 1:size(spirals_raw)
    frame = spirals_raw(i,5);
    spiralx = spirals_raw(i,1); spiraly = spirals_raw(i,2);

    spiral_index = find(spirals_predict(:,5)==frame);
    spiral_p1 = spirals_predict(spiral_index,:);
    index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
        spiral_p1(:,1)<= spiralx+100 & ...
        spiral_p1(:,2)>= spiraly-100 & ...
        spiral_p1(:,2)<= spiraly+100);
    spiral_p2 = spiral_p1(index1,:);
    if size(spiral_p2,1) == 1
        spiral_temp = [spirals_raw(i,1:5), spiral_p2];
        spiral_match = [spiral_match;spiral_temp];
    end
end
spiral_match(:,11) = 0;
match = spiral_match(:,9)-spiral_match(:,4);                 % check if spiral directions are matching in raw and predicted
match = not(match);
spiral_match(:,11) = match;