function getSpiralsCountsBySession(T_all,spiral_all,radius)
labels = {"correct","incorrect","miss"};
%% 
for kk = 1:3
%% left stim
    clear index spiral_cell_temp spiral_cell_temp2
    index = (T_all.label == labels{kk} & T_all.left_contrast > 0);
    spiral_cell_temp = spiral_all(index,:);
    trialN1 = sum(index);
    spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
    for n = 1:141
        clear spiral_temp indx1 indx2
        spiral_temp =  spiral_cell_temp2{n};
        indx1 = (spiral_temp(:,1) > 570 & spiral_temp(:,3)==100);
        spiral_Lstim_contra(n,1) = sum(indx1);
        indx2 = (spiral_temp(:,1) < 570 & spiral_temp(:,3)==100);
        spiral_Lstim_ipsi(n,1) = sum(indx2);
    end
    %% right stim
    clear index spiral_cell_temp spiral_cell_temp2
    index = (T_all.label == labels{kk} & T_all.right_contrast > 0);
    spiral_cell_temp = spiral_all(index ,:);
    trialN2 = sum(index);
    spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
    for n = 1:141
        clear spiral_temp indx1 indx2
        spiral_temp =  spiral_cell_temp2{n};
        indx1 = (spiral_temp(:,1) < 570 & spiral_temp(:,3)==100);
        spiral_Rstim_contra(n,1) = sum(indx1);
        indx2 = (spiral_temp(:,1) > 570 & spiral_temp(:,3)==100);
        spiral_Rstim_ipsi(n,1) = sum(indx2);
    end
    %%
    spiral_contra(:,i,kk) = (spiral_Lstim_contra+spiral_Rstim_contra)./(trialN1+trialN2);
    spiral_ipsi(:,i,kk) = (spiral_Lstim_ipsi+spiral_Rstim_ipsi)./(trialN1+trialN2);
end
%%
clear index spiral_cell_temp spiral_cell_temp2
index = (T_all.label == "miss");
spiral_cell_temp = spiral_all(index,:);
trialN = sum(index);
spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
for n = 1:141
    clear spiral_temp indx1 indx2
    spiral_temp =  spiral_cell_temp2{n};
    indx1 = (spiral_temp(:,3)==100);
    spiral_count_sum_miss(n,i) = sum(indx1)/trialN*1;
end    