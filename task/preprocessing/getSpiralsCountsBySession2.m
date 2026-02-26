function [spiral_count_pre, spiral_count_post] = getSpiralsCountsBySession2(T_all,spiral_all,radius)
%%
labels = {"correct","incorrect","miss"};
spiral_count_sum_all =  zeros(3,141);
for id = 1:3
    clear index sprial_label
    label = labels{id};   
    index = (T_all.label == label & abs(T_all.left_contrast-T_all.right_contrast)> 0);
    spiral_label = spiral_all(index ,:);
    trialN = sum(index);
    %%
    spiral_count_sum =zeros(141,1);
    for n = 1:141
        clear spiral_temp spiral_temp2 indx 
        spiral_temp = spiral_label(:,n);
        spiral_temp2 = cat(1,spiral_temp{:});
        if not(isempty(spiral_temp2))
            indx = ismember(spiral_temp2(:,3),radius);
        else 
            indx = 0;
        end
        spiral_count_sum(n,1) = sum(indx)/trialN*1;
    end
    spiral_count_sum_all(id,:) = spiral_count_sum;
end
%%
spiral_count_sum_all = spiral_count_sum_all*35;
% spiral_count_pre = squeeze(sum(spiral_count_sum_all(:,63:70),2))/8;
% spiral_count_post = squeeze(sum(spiral_count_sum_all(:,75:82),2))/8;

% spiral_count_pre = squeeze(sum(spiral_count_sum_all(:,63:67),2))/5;
% spiral_count_post = squeeze(sum(spiral_count_sum_all(:,77:81),2))/5;

spiral_count_pre = squeeze(sum(spiral_count_sum_all(:,64:66),2))/3;
spiral_count_post = squeeze(sum(spiral_count_sum_all(:,78:80),2))/3;