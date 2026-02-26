function [predicted_signals,explained_var] = ...
    kernelPrediction(k_cv,regressors,signals,t_shifts,lambdas,zs,cvfold,return_constant,use_constant,discontinuities)

% Use no temporal delay if none specified
if ~exist('t_shifts','var') || isempty(t_shifts)
    t_shifts = 0;
end

% Convert regressors and t_shifts to cells if not already
if ~iscell(regressors)
    regressors = {regressors};
end
if ~iscell(t_shifts)
    t_shifts = {t_shifts};
end

% Standardize orientations
regressors = reshape(regressors,1,[]);
t_shifts = reshape(t_shifts,1,[]);
lambdas = reshape(lambdas,1,[]);

% Z-score all regressors and signals to get beta weights if selected
if ~exist('zs','var') || isempty(zs)
    zs = [false,true];
end
if zs(1)
    regressors = cellfun(@(x) zscore(x,[],2),regressors,'uni',false);
end
if zs(2)
    signals = zscore(signals,[],2);
end

% Set cross-validation to 1-fold if not entered
if ~exist('cvfold','var') || isempty(cvfold)
   cvfold = 1; 
end

% Set return_constant false if not entered
if ~exist('return_constant','var') || isempty(return_constant)
   return_constant = false; 
end

% Set use_constant true if not entered
if ~exist('use_constant','var') || isempty(use_constant)
   use_constant = true; 
end

% Set/check discontinuities
if ~exist('discontinuities','var') || isempty(discontinuities)
   discontinuities = zeros(1,size(signals,2)); 
elseif length(discontinuities) ~= size(signals,2)
    error('Discontinuities vector doesn''t match signals length')        
end
% (the binning and end are always discontinuities)
discontinuities([1,end]) = 1;
% (force discontinuities to be a column matrix)
discontinuities = reshape(discontinuities,[],1);

% Create design matrix of all time-shifted regressors
regressor_design = cellfun(@(regressors,t_shifts) repmat(regressors', ...
    [1,1,length(t_shifts)]),regressors,t_shifts,'uni',false);

% Temporally shift each page (regressors and discontinuities)
for curr_regressors = 1:length(regressor_design)
    
    % Set up discontinuities matrix (none at t shift == 0)
    curr_discontinuities = repmat(discontinuities,1, ...
        size(regressor_design{curr_regressors},2), ...
        size(regressor_design{curr_regressors},3));
    curr_discontinuities(:,:,t_shifts{curr_regressors} == 0) = 0;
    
    for curr_kernel_frame = 1:length(t_shifts{curr_regressors})
        regressor_design{curr_regressors}(:,:,curr_kernel_frame) = ...
            circshift(regressor_design{curr_regressors}(:,:,curr_kernel_frame), ...
            [t_shifts{curr_regressors}(curr_kernel_frame),0,0]);
        
        curr_discontinuities(:,:,curr_kernel_frame) = ...
            circshift(curr_discontinuities(:,:,curr_kernel_frame), ...
            [t_shifts{curr_regressors}(curr_kernel_frame),0,0]);
    end
    
    % Cumulative sum the discontinuities in appropriate direction
    curr_discontinuities_cumulative = curr_discontinuities > 0;
    curr_discontinuities_cumulative(:,:,t_shifts{curr_regressors} < 0) = ...
        cumsum(curr_discontinuities(:,:,t_shifts{curr_regressors} < 0),3,'reverse') > 0;
    curr_discontinuities_cumulative(:,:,t_shifts{curr_regressors} > 0) = ...
        cumsum(curr_discontinuities(:,:,t_shifts{curr_regressors} > 0),3) > 0;
    
    % Zero regressors predicting discontinuous / invalid locations
    regressor_design{curr_regressors}(curr_discontinuities_cumulative) = 0;
    
end

regressor_design = cell2mat(cellfun(@(regressor_design) ...
    reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)), ...
    regressor_design,'uni',false));

% Ridge regression for reducing noise: add offsets to design matrix to penalize k
if exist('lambdas','var') && any(lambdas)
    if length(lambdas) == 1
        ridge_matrix = lambdas*eye(size(regressor_design,2));
    elseif length(lambdas) == length(regressors)
        lambda_vector = cell2mat(reshape(cellfun(@(reg,t,lam) repmat(lam,size(reg,1)*length(t),1), ...
            regressors,t_shifts,num2cell(lambdas),'uni',false),[],1));
        ridge_matrix = bsxfun(@times,eye(size(regressor_design,2)),lambda_vector);
    else
        error('Number of lambdas doesn''t match regressor groups');
    end
else
    ridge_matrix = [];
end

% Prepare column of 1's to have a constant term (if selected)
if use_constant
    constant = ones(size(regressor_design,1),1);
    % if there's a ridge matrix, add another row and column of zeros
    if ~isempty(ridge_matrix)
        ridge_matrix(end+1,end+1) = 0;
    end
else
    constant = [];
end

regressors_gpu = gpuArray([[regressor_design,constant];ridge_matrix]);

% Regression (and cross validation if selected)

% Check for NaNs: 
% Any signals which are all NaN won't be regressed
predictable_signals = ~all(isnan(signals),2);
% If there are NaNs not shared across the remaining signals, error out
% if ~all(all(isnan(signals(predictable_signals,:)),1) | ...
%         all(~isnan(signals(predictable_signals,:)),1))
%     error('NaN values vary across signals')
% end
    
% Predictable samples: non-zero/nan regressors, non-nan signals
predictable_samples = ...
    any(regressor_design,2) & ...
    ~any(isnan(regressor_design),2) & ...
    ~any(isnan(signals(predictable_signals,:)),1)';

predicted_signals = nan(size(signals));
for curr_cv = 1:cvfold
    
    train_idx = predictable_samples;
    test_idx = predictable_samples;
   
    % If regressors for training fold are empty, warning
    if ~all(any(regressors_gpu(train_idx,:),1))
        warning('Regressors in fold unfilled (not enough trials?)');
    end

    predicted_signals(predictable_signals,test_idx) = ...
        gather(regressors_gpu(test_idx,:)*gpuArray(k_cv(:,predictable_signals,curr_cv)))';
end

% Throw errors for bad numbers in kernel or predicted signals
if ...
        any(reshape(isnan(k_cv(:,predictable_signals,:)),[],1)) || ...
        any(reshape(isnan(predicted_signals(predictable_signals,predictable_samples)),[],1)) || ...
        any(reshape(isinf(k_cv(:,predictable_signals,:)),[],1)) || ...
        any(reshape(isinf(predicted_signals(predictable_signals,:)),[],1)) ...
    error('Inf/NaN in kernel or predicted signal')
end

% Total explained variance (R^2)
sse_residual = sum((signals(predictable_signals,predictable_samples) - ...
    predicted_signals(predictable_signals,predictable_samples)).^2,2);
sse_total = sum((signals(predictable_signals,predictable_samples) - ...
    nanmean(signals(predictable_signals,predictable_samples),2)).^2,2);
explained_var.total = nan(size(signals,1),1);
explained_var.total(predictable_signals) = 1 - (sse_residual./sse_total);









