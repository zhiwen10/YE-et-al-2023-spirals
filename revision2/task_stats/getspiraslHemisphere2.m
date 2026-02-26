function spiral = getspiraslHemisphere2(T_session,spiral_session,label)
%% left stim
spiral_Lstim_contra = [];
spiral_Lstim_ipsi = [];
clear index spiral_cell_temp spiral_cell_temp2
index1 = (T_session.label == label & T_session.left_contrast > 0); %left stim
index2 = (T_session.label == label & T_session.right_contrast > 0); % right_stim
index_all = {index1;index2};
% left stim is default, flip spirals for right stim
for k  = 1:2
    index = index_all{k};
    spiral_cell_temp = spiral_session(index,:);
    if not(isempty(spiral_cell_temp))
        trialN1 = sum(index);
        spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
        %%
        for n = 1:141
            clear spiral_temp spiral_temp2a spiral_temp2b indx1 indx2
            spiral_temp =  spiral_cell_temp2{n};
            if k == 2
                spiral_temp(:,1) = 1140-spiral_temp(:,1);
                spiral_temp(:,4) = -spiral_temp(:,4);
            end
            %%
            % for left stim, right hemi is contralateral side
            indx1 = (spiral_temp(:,1) > 570 );
            spiral_temp2a = spiral_temp(indx1,:);
            if not(isempty(indx1))
                spiral_Lstim_contra = [spiral_Lstim_contra;spiral_temp2a];
            end
            % for left stim, left hemi is ipsilateral side
            indx2 = (spiral_temp(:,1) < 570);
            spiral_temp2b = spiral_temp(indx2,:);
            if not(isempty(indx2))
                spiral_Lstim_ipsi = [spiral_Lstim_ipsi;spiral_temp2b];
            end
        end
    end
end
%%
[lia1,locb1] = ismember(spiral_Lstim_contra(:,1:2),indexSSp_right,'rows');
spirals_post_temp2 = spiral_Lstim_contra(lia1,:);
ratio_contra = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
%%
[lia1,locb1] = ismember(spiral_Lstim_ipsi(:,1:2),indexSSp_left,'rows');
spirals_post_temp2 = spiral_Lstim_ipsi(lia1,:);
ratio_ipsi = sum(spirals_post_temp2(:,4) == -1)./size(spirals_post_temp2,1);
