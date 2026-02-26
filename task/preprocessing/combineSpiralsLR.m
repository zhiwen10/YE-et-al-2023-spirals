function spiral_all = combineSpiralsLR(spiral_L, spiral_R)
spiral_all = cell(141,1);
    for i = 1:141
        clear spiral_L_temp spiral_R_temp spiral_temp ...
        %% combine correct left and right
        spiral_L_temp = spiral_L{i};
        spiral_R_temp = spiral_R{i};
        spiral_temp = cat(1,spiral_L_temp,spiral_R_temp); 
        spiral_all{i,1} = spiral_temp;    
    end