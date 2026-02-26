function [ratio_right,ratio_left] = getspiraslHemisphere3(T_session,spiral_session,label,radius,indexSSp_right,indexSSp_left)
%% left stim
spiral = [];
spiral_Lstim_ipsi = [];
clear index spiral_cell_temp spiral_cell_temp2
index1 = (T_session.label == label & T_session.left_contrast > 0); %left stim
index2 = (T_session.label == label & T_session.right_contrast > 0); % right_stim
index_all = {index1;index2};
%%
% left stim is default, flip spirals for right stim
spiral = cell(141,2);
for k  = 1:2
    %%
    index = index_all{k};
    spiral_cell_temp = spiral_session(index,:);
    %%
    if not(isempty(spiral_cell_temp))
        trialN1 = sum(index);
        spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
        %%
        for n = 1:141
            clear spiral_temp
            spiral_temp =  spiral_cell_temp2{n};
            if not(isempty(spiral_temp))
                if k == 2
                    spiral_temp(:,1) = 1140-spiral_temp(:,1);
                    spiral_temp(:,4) = -spiral_temp(:,4);
                end           
                spiral{n,k} = spiral_temp;
            end
        end
    end
end
%%
if not(isempty(spiral))
    spiral2 = cell(141,1);
    for i = 1:141
        clear cspiral_temp2
        spiral_temp2 = spiral(i,:);
        spiral2{i,1} = cat(1,spiral_temp2{:});
    end
    frames_pre = 63:67; frames_post = 77:81;
    spirals3a = cat(1,spiral2{frames_pre});
    spirals3b = cat(1,spiral2{frames_post});
    spirals4 = {spirals3a,spirals3b};
    %%
    for i = 1:2 % pre and post
        spiral_current = spirals4{i};
        if not(isempty(spiral_current))
            if not(isempty(radius))
                spiral_current = spiral_current(spiral_current(:,3)== radius,:);
            end
            [lia1,locb1] = ismember(spiral_current(:,1:2),indexSSp_right,'rows');
            spirals_post_temp2 = spiral_current(lia1,:);
            ratio_right(i,1) = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
            [lia1,locb1] = ismember(spiral_current(:,1:2),indexSSp_left,'rows');
            spirals_post_temp2 = spiral_current(lia1,:);
            ratio_left(i,1) = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
        else
            ratio_right(i,1) = nan;
            ratio_left(i,1) = nan;
        end
    end
else
    ratio_right = nan(2,1);
    ratio_left = nan(2,1);
end