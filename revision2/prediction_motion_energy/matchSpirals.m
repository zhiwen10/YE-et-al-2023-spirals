function [spiral_left_match,spiral_right_match] =  matchSpirals(spirals_filt,brain_index_left,brain_index_right,pwAll,tform,traceAmp)
%% search for match in left and right hemispehre seperately
[liaL,locbL] = ismember(spirals_filt(:,1:2),brain_index_left,'rows');
[liaR,locbR] = ismember(spirals_filt(:,1:2),brain_index_right,'rows');
spirals_left = spirals_filt(liaL,:);
spirals_right = spirals_filt(liaR,:);
%% load spirals from prediction
[spiralsT1(:,1),spiralsT1(:,2)] = transformPointsForward(...
    tform,pwAll(:,1),pwAll(:,2));
pwAll(:,1:2) = round(spiralsT1); 
[liaL,locbL] = ismember(pwAll(:,1:2),brain_index_left,'rows');
spirals_prediction_left = pwAll(liaL,:);
[liaR,locbR] = ismember(pwAll(:,1:2),brain_index_right,'rows');
spirals_prediction_right = pwAll(liaR,:);
%% find spirals within 100 pixels distance from each raw spiral
% only care about the situation when only one spiral is within this
% distance
clear match spiral_p1 spiral_p2 index1 sprial_index
spiral_left_match = [];                                                % check left hemisphere first 
for i = 1:size(spirals_left)
    frame = spirals_left(i,5);
    spiralx = spirals_left(i,1); spiraly = spirals_left(i,2);

    spiral_index = find(spirals_prediction_left(:,5)==frame);
    spiral_p1 = spirals_prediction_left(spiral_index,:);
    index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
        spiral_p1(:,1)<= spiralx+100 & ...
        spiral_p1(:,2)>= spiraly-100 & ...
        spiral_p1(:,2)<= spiraly+100);
    spiral_p2 = spiral_p1(index1,:);
    if size(spiral_p2,1) == 1
        spiral_temp = [spirals_left(i,1:5), spiral_p2];
        spiral_left_match = [spiral_left_match;spiral_temp];
    end
end
spiral_left_match(:,11) = 0;
match = spiral_left_match(:,9)-spiral_left_match(:,4);                 % check if spiral directions are matching in raw and predicted
match = not(match);
spiral_left_match(:,11) = match;
spiral_left_match(:,12) = traceAmp(1,spiral_left_match(:,5));          % attach 2-8Hz amplitude value for that spiral frame
%%
clear match spiral_p1 spiral_p2 index1 sprial_index
spiral_right_match = [];                                               % now check right hemisphere 
for i = 1:size(spirals_right)
    frame = spirals_right(i,5);
    spiralx = spirals_right(i,1); spiraly = spirals_right(i,2);

    spiral_index = find(spirals_prediction_right(:,5)==frame);
    spiral_p1 = spirals_prediction_right(spiral_index,:);
    index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
        spiral_p1(:,1)<= spiralx+100 & ...
        spiral_p1(:,2)>= spiraly-100 & ...
        spiral_p1(:,2)<= spiraly+100);
    spiral_p2 = spiral_p1(index1,:);
    if size(spiral_p2,1) == 1
        spiral_temp = [spirals_right(i,1:5), spiral_p2];
        spiral_right_match = [spiral_right_match;spiral_temp];
    end
end
spiral_right_match(:,11) = 0;
match = spiral_right_match(:,9)-spiral_right_match(:,4);               % check if spiral directions are matching in raw and predicted
match = not(match);
spiral_right_match(:,11) = match;
spiral_right_match(:,12) = traceAmp(1,spiral_right_match(:,5));        % attach 2-8Hz amplitude value for that spiral frame