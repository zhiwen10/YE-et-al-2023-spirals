function [T_ratio,rt_median,trial_count,reaction_time_sort] = sort_ratio_by_contrast2(T_all)
contrast = [0.06,0.125,0.25,0.5,1]'; 
%%
ratio = [];
reaction_time_sort = cell(5,1);
T_sort = {};
 for i = 1:5
     clear T_temp
     T_temp = T_all(abs(T_all.left_contrast- T_all.right_contrast)== contrast(i),:);
     ratio(i,:) = sum(T_temp{:,5:10},1)./size(T_temp,1);
     T_sort{i} = T_temp;
     reaction_time_sort{i} = T_temp(ismember(T_temp.label,{'correct','falarmL','falarmR'}),:).wheel_onset;
     trial_count(i,1) = size(T_temp,1);
 end
%%
rt_median = cellfun(@mean, reaction_time_sort)'; 
%%
ratio = [contrast,ratio];
T_ratio = array2table(ratio);
T_ratio = renamevars(T_ratio,["ratio1","ratio2","ratio3","ratio4","ratio5","ratio6","ratio7"], ...
                 ["contrast","miss","correct","incorrect","reject","falarmL","falarmR"]);