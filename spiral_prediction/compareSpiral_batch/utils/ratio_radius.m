
function [mean_ratio,std_ratio,sem_ratio,ratio_all1] = ratio_radius(T1,spiral_left_match_all,spiral_right_match_all)
for kk = 1:size(T1,1)
    spiral_left_temp = spiral_left_match_all{kk};
    spiral_right_temp = spiral_right_match_all{kk};
    count = 1;
    for radius = 40:10:100
        match = []; spiral_all_match_temp = [];
        indx1 = (spiral_left_temp(:,3) == radius);
        spiral_left_match_temp = spiral_left_temp(indx1,:);
        indx2 = (spiral_right_temp(:,3) == radius);
        spiral_right_match_temp = spiral_right_temp(indx2,:);
        spiral_all_match_temp = [spiral_left_match_temp;spiral_right_match_temp];
        if size(spiral_all_match_temp,1)<100
            ratio_all1(count,kk) = nan;
        else         
            match = spiral_all_match_temp(:,9)-spiral_all_match_temp(:,4);
            match = not(match);
            ratio1 = sum(match)/size(spiral_all_match_temp,1);
            ratio_all1(count,kk) = ratio1;
        end
        count = count+1;
    end
end
mean_ratio= mean(ratio_all1,2,"omitnan");
std_ratio = std(ratio_all1,[],2,"omitnan");
count_all = sum(not(isnan(ratio_all1)),2);
sem_ratio = std_ratio./sqrt(count_all);