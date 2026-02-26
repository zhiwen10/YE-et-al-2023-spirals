function [dV_predict,k_cv] = mua_prediction_full(dV1_train,MUA_std_train,dV1_test,MUA_std_test,kernel_t,varargin)
%%
if nargin == 4
    perm = varargin{1};
else
    perm = 0;
end
% Set parameters for regression
use_svs = 1:50;
sample_rate = 35;
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [true,false];
cvfold = 5;
%% train on odd epochs
lambda1 = 1;
return_constant = 1;
for i = 1:numel(lambda1)
    lambda = lambda1(i);
    [k{i},predicted_spikes(:,:,i),explained_var{i}] = ...
        AP_regresskernel_mod(MUA_std_train,dV1_train, ...
        kernel_frames,lambda,zs,cvfold,return_constant);
end
%%
k_cv = k{1}{1};
k_cv = reshape(k_cv,size(k_cv,1)*size(k_cv,2),size(k_cv,3));
k_cv(end+1,:) = squeeze(k{1}{2});
%%
cvfold = 1;
for i = 1:5
    MUA_std_test1 = squeeze(MUA_std_test(:,:,i));
    dV1_test1 = squeeze(dV1_test(:,:,i));
    %% 
    [dV_predict(:,:,i),explained_var2(:,i)] = ...
        kernelPrediction(k_cv,MUA_std_test1,dV1_test1,kernel_frames,lambda1,zs,cvfold);
end      
end            
