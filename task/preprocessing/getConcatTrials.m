function [spiral_all2] = getConcatTrials(spiral_all)
spiral_all2 = cell(141,1);
for i = 1:141
    clear current_correct
    current_correct = spiral_all(:,i);
    spiral_all2{i,1} = cat(1,current_correct{:});
end