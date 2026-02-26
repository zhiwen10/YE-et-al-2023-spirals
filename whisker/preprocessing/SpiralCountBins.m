function spiral_count_sum_all = SpiralCountBins(spiral_cell)
%%
radius = 40:10:100;
trialN = size(spiral_cell,1);
for j = 1:7
    spiral_count =zeros(size(spiral_cell));
    for m = 1:size(spiral_cell,1)
        spiral_cell_temp = spiral_cell(m,:);
        %%
        for n = 1:size(spiral_cell,2)
            spiral_temp = spiral_cell_temp{n};
            if not(isempty(spiral_temp))
                indx = (spiral_temp(:,3)== radius(j));
                spiral_temp2 = spiral_temp(indx,:);
                if not(isempty(spiral_temp2))
                    spiral_count(m,n) = 1;
                end            
            end
        end
    end
    spiral_count_sum = sum(spiral_count,1)./trialN;
    spiral_count_sum_all(j,:) = spiral_count_sum;
end