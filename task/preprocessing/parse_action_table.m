function [T,T_ratio] = parse_action_table(block)    
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);                 % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1);
% get response no repeat
response = response(norepeat_indx)';
% get feedback no repeat
feedback = block.events.feedbackValues(norepeat_indx)';
% get stim no repeat
left_contrast = block.events.contrastLeftValues;
right_contrast = block.events.contrastRightValues;
left_contrast = left_contrast(norepeat_indx)';
right_contrast = right_contrast(norepeat_indx)';
% reactiom time
reaction_time = block.events.responseTimes(1:ntrial)-block.events.stimulusOnTimes(1:ntrial);
reaction_time = reaction_time(norepeat_indx)';
%% get table
T = table(left_contrast,right_contrast,response,feedback);
% miss
miss = zeros(size(T,1),1);
miss(logical(T.right_contrast-T.left_contrast) & T.response == 0,1)= 1;
T.miss = miss;
% correct
correct = zeros(size(T,1),1);
correct((logical(T.right_contrast-T.left_contrast)& T.feedback == 1),1)= 1;
T.correct = correct;
% incorrect
incorrect = zeros(size(T,1),1);
% incorrect((logical(contrastRight-contrastLeft)& not(miss) & (choice~=action)),1)= 1;
incorrect((logical(T.right_contrast-T.left_contrast)& logical(T.response) & T.feedback==0),1)= 1;
T.incorrect = incorrect;
% reject
reject = zeros(size(T,1),1);
reject((T.right_contrast == 0 & T.left_contrast == 0 & T.response==0),1)= 1;
T.reject = reject;
% falarm
falarmL = zeros(size(T,1),1);
falarmL((T.right_contrast == 0 & T.left_contrast == 0& T.response == 1),1) = 1;
T.falarmL = falarmL;

falarmR = zeros(size(T,1),1);
falarmR((T.right_contrast == 0 & T.left_contrast == 0& T.response == -1),1) = 1;
T.falarmR = falarmR;
% check if all trials were conted

sum_T = sum(T{:,5:10},2);
match = all(sum_T);
T_label = T(:,5:10);
label1 = {};
list = {'miss','correct','incorrect','reject','falarmL','falarmR'};
for i = 1:size(T_label,1)
    binary_temp = T_label{i,:};
    a = list{find(binary_temp)};
    label1(i,1) = {a};
end
T.label = label1;
T.reactionTime = reaction_time;
%%
T.left_contrast(T.left_contrast==0.006) = 0.06;
T.right_contrast(T.right_contrast==0.006) = 0.06;
%%
contrast = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1]'; 
stim = [1, 0; 0.5, 0; 0.25, 0;0.125, 0; 0.06,0; 0, 0; 0,0.06; 0,0.125; 0, 0.25; 0, 0.5; 0, 1];
ratio = [];
reaction_time_sort = {};
T_sort = {};
 for i = 1:11
     clear T_temp
     T_temp = T(T.left_contrast==stim(i,1) & T.right_contrast==stim(i,2),:);
     ratio(i,:) = sum(T_temp{:,5:10},1)./size(T_temp,1);
     T_sort{i} = T_temp;
     reaction_time_sort{i} = T_temp(ismember(T_temp.label,{'correct','falarmL','falarmR'}),:).reactionTime;
 end
%%
for i = 1:11
    rt_temp = reaction_time_sort{i};
    rt_mean(i,1) = mean(rt_temp);
    rt_sem(i,1) = std(rt_temp)./sqrt(size(rt_temp,1));
end
%%
ratio = [contrast,ratio];
T_ratio = array2table(ratio);
T_ratio = renamevars(T_ratio,["ratio1","ratio2","ratio3","ratio4","ratio5","ratio6","ratio7"], ...
                 ["contrast","miss","correct","incorrect","reject","falarmL","falarmR"]);