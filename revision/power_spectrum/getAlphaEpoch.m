function [alpha_epoch2,alpha_binary] = getAlphaEpoch(traceAmp2,threshold)
%%
% threshold = 0.0025;
alpha1 = (traceAmp2>threshold);
alpha1_diff = diff(alpha1);
alpha1_diff = [0,alpha1_diff];
% plot(t, alpha2_diff/100,'c');
epoch_start = find(alpha1_diff==1);
epoch_end = find(alpha1_diff==-1);

if epoch_start(1)<epoch_end(1)
    epoch_n = min(numel(epoch_start),numel(epoch_end));
    alpha_epoch = [epoch_start(1:epoch_n)',epoch_end(1:epoch_n)']; 
    epoch_length = alpha_epoch(:,2)-alpha_epoch(:,1);
    alpha_epoch2 = alpha_epoch(epoch_length>17,:);
    alpha_binary = zeros(size(alpha1));
else
    epoch_start = [1,epoch_start];
    epoch_n = min(numel(epoch_start),numel(epoch_end));
    alpha_epoch = [epoch_start(1:epoch_n)',epoch_end(1:epoch_n)']; 
    epoch_length = alpha_epoch(:,2)-alpha_epoch(:,1);
    alpha_epoch2 = alpha_epoch(epoch_length>17,:);
    alpha_binary = zeros(size(alpha1));
end
for i = 1:size(alpha_epoch2,1)
    alpha_binary(alpha_epoch2(i,1):alpha_epoch2(i,2)) = 1;
end
alpha_binary = logical(alpha_binary);
% alpha_ratio = sum(alpha_binary)./numel(alpha_binary);
%%
