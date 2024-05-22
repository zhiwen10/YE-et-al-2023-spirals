function [dV_raw,dV_predict,epoch_indx] = get_prediction(dV1,MUA_std,varargin)
%%
if nargin == 3
    perm = varargin{1};
else
    perm = 0;
end
% Set parameters for regression
use_svs = 1:50;
% kernel_t = [-3/35, 3/35];
kernel_t = [-0.2,0.2];
sample_rate = 35;
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [true,false];
cvfold = 5;
%% get 10 segment sample index
sample_length = size(dV1,2);
epochN = 10;
epochSize = floor(sample_length/epochN);
for i = 1:epochN
    epoch_indx(i,1) = 1+epochSize*(i-1);
    epoch_indx(i,2) = epochSize*i;
end
interval = epoch_indx(:,2)-epoch_indx(:,1);
%% train on odd epochs
% lambda1 = (1:10:100).^2;
lambda1 = 1;
return_constant = 1;
train_indx = [];
for i = 1:5
    train_indx1 = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
    train_indx = [train_indx,train_indx1];
end
for i = 1:numel(lambda1)
    lambda = lambda1(i);
    [k{i},predicted_spikes(:,:,i),explained_var{i}] = ...
        AP_regresskernel_mod(MUA_std(:,train_indx),dV1(:,train_indx), ...
        kernel_frames,lambda,zs,cvfold,return_constant);
end
%% predict even epochs seperately
clear test_indx

% lambda1 = 121;
cvfold = 1;
for i = 1:5
    test_indx = [epoch_indx(2*i,1):epoch_indx(2*i,2)];
    k_cv = k{1}{1};
    k_cv = reshape(k_cv,size(k_cv,1)*size(k_cv,2),size(k_cv,3));
    k_cv(end+1,:) = squeeze(k{1}{2});
%%   
if perm
    randIndx = randperm(size(MUA_std,1));
    MUA_std_test = MUA_std(randIndx,test_indx);
else
    MUA_std_test = MUA_std(:,test_indx);
end
%% 
    [dV_predict(:,:,i),explained_var2(:,i)] = ...
        kernelPrediction(k_cv,MUA_std_test,dV1(:,test_indx),kernel_frames,lambda1,zs,cvfold);
    dV_raw(:,:,i) = dV1(:,test_indx);
end            
