function [epoch_indx,train_indx,test_indx] = divide_epoch(dV1, epochN)
% dV1 = dimensions  x time_series
% divide by second dimension
sample_length = size(dV1,2);
% epochN = 10;
epochSize = floor(sample_length/epochN);
for i = 1:epochN
    epoch_indx(i,1) = 1+epochSize*(i-1);
    epoch_indx(i,2) = epochSize*i;
end

% concatenate all odd epochs
train_indx = [];
for i = 1:5
    train_indx1 = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
    train_indx = [train_indx,train_indx1];
end

test_indx = epoch_indx(2:2:end,:);
end