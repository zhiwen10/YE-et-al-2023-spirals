function [N,edges,phase_mu,phase_var,flow_mu,flow_var,traceAmp_mean] = compare_flow1(Ut,mimg1,dV1,MUA_std,BW,len)
%%
% len is size of sample to use 
%%
% Set parameters for regression
params.use_svs = 1:50;
kernel_t = [-0.5,0.5];
sample_rate = 35;
params.kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
params.zs = [true,false];
params.cvfold = 5;
params.lambda = 1;
params.return_constant = 1;
flow =1; 
roi = 1;
%% get 10 segment sample index
epochN = 10;
[epoch_indx,train_indx,test_indx] = divide_epoch(dV1, epochN);
interval = epoch_indx(:,2)-epoch_indx(:,1);
%% train on odd epochs
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel_mod(MUA_std(:,train_indx),dV1(:,train_indx), ...
    params.kernel_frames,params.lambda,params.zs,params.cvfold,params.return_constant);
%% predict even epochs seperately
cvfold = 1;
% use 2nd batch in test batch
ibatch = 2;
test_indx_i = [test_indx(ibatch,1):test_indx(ibatch,2)];
k_cv = k{1};
k_cv = reshape(k_cv,size(k_cv,1)*size(k_cv,2),size(k_cv,3));
k_cv(end+1,:) = squeeze(k{2});
%% raw phase and flow
dV_raw = dV1(:,test_indx_i);  
[tracePhase_raw,vxy_raw,traceAmp_raw] =get_flowfield5(Ut,dV_raw,mimg1,flow,len);
vxy_raw(:,~BW) = nan;
traceAmp_raw(:,~BW) = nan;
tracePhase_raw(:,~BW) = nan;
traceAmp_mean = mean(traceAmp_raw,[2,3],"omitnan");
%% predicted phase and flow and scramble
MUA_std_test = MUA_std(:,test_indx_i);
dV1_test = dV1(:,test_indx_i);
for count = 1:21
    %% from count2 start to scramble
    if count>1
        randIndx = randperm(size(MUA_std,1));
        MUA_std_test = MUA_std_test(randIndx,:);
    end
    
    [dV_predict,explained_var2] = kernelPrediction(k_cv,MUA_std_test,dV1_test,params.kernel_frames,params.lambda,params.zs,cvfold);
    [tracePhase_predict,vxy_predict,traceAmp_predict] =get_flowfield5(Ut,dV_predict,mimg1,flow,len);   
    vxy_predict(:,~BW) = nan;
    tracePhase_predict(:,~BW) = nan;
    %%
    [N(count,:),edges(count,:),phase_mu(count,:),phase_var(count,:),flow_mu(count,:),flow_var(count,:)] = get_flow_metric1(traceAmp_raw,tracePhase_raw,tracePhase_predict,vxy_raw,vxy_predict);
end
end            