%%
function results = permutation_test_paired_timepoints(data, varargin)
    % PERMUTATION_TEST_PAIRED_TIMEPOINTS Permutation test for paired time point comparison
    %
    % Inputs:
    %   data - struct with fields:
    %          .Subject - vector of subject IDs
    %          .TimePoint - vector of time points (1 or 2)
    %          .Value - vector of measured values
    %          .Session - vector of session IDs (optional, for tracking)
    %
    % Optional parameters (name-value pairs):
    %   'nPermutations' - number of permutations (default: 10000)
    %   'randomSeed' - random seed for reproducibility (default: 42)
    %   'showProgress' - show progress updates (default: true)
    %   'makePlots' - create visualization plots (default: true)
    %
    % Returns:
    %   results - struct with test results
    
    % Parse input arguments
    p = inputParser;
    addRequired(p, 'data', @isstruct);
    addParameter(p, 'nPermutations', 10000, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'randomSeed', 42, @isnumeric);
    addParameter(p, 'showProgress', true, @islogical);
    addParameter(p, 'makePlots', true, @islogical);
    parse(p, data, varargin{:});
    
    nPermutations = p.Results.nPermutations;
    randomSeed = p.Results.randomSeed;
    showProgress = p.Results.showProgress;
    makePlots = p.Results.makePlots;
    
    % Set random seed for reproducibility
    rng(randomSeed);
    
    % Validate input data
    validateInputData(data);
    
    % Calculate observed test statistic
    tp1_indices = data.TimePoint == 1;
    tp2_indices = data.TimePoint == 2;
    tp1_values = data.Value(tp1_indices);
    tp2_values = data.Value(tp2_indices);
    
    % observed_effect = mean(tp2_values) - mean(tp1_values);
    observed_effect = mean(tp2_values - tp1_values);
    
    if showProgress
        fprintf('Observed effect (TP2 - TP1): %.4f\n', observed_effect);
        fprintf('Sample sizes: TP1=%d, TP2=%d\n', length(tp1_values), length(tp2_values));
    end
    
    % Get unique subjects
    subjects = unique(data.Subject);
    n_subjects = length(subjects);
    
    if showProgress
        fprintf('Number of subjects: %d\n', n_subjects);
        fprintf('Running %d permutations...\n', nPermutations);
    end
    
    % Preallocate for efficiency
    permuted_effects = zeros(nPermutations, 1);
    
    % Main permutation loop
    for perm = 1:nPermutations
        % Initialize for this permutation
        perm_tp1_values = [];
        perm_tp2_values = [];
        
        % For each subject, randomly reassign timepoint labels
        for s = 1:n_subjects
            subject_idx = data.Subject == subjects(s);
            subject_values = data.Value(subject_idx);
            subject_timepoints = data.TimePoint(subject_idx);
            
            % Get original session counts per timepoint
            tp1_count = sum(subject_timepoints == 1);
            tp2_count = sum(subject_timepoints == 2);
            total_sessions = tp1_count + tp2_count;
            
            % Randomly select which sessions go to TP1
            tp1_indices = randperm(total_sessions, tp1_count);
            tp2_indices = setdiff(1:total_sessions, tp1_indices);
            
            % Assign values to permuted timepoints
            perm_tp1_values = [perm_tp1_values; subject_values(tp1_indices)];
            perm_tp2_values = [perm_tp2_values; subject_values(tp2_indices)];
        end
        
        % Calculate permuted effect
        % permuted_effects(perm) = mean(perm_tp2_values) - mean(perm_tp1_values);
        permuted_effects(perm) = mean(perm_tp2_values - perm_tp1_values);
        
        % Progress indicator
        if showProgress && mod(perm, 1000) == 0
            fprintf('Completed %d/%d permutations\n', perm, nPermutations);
        end
    end
    
    % Calculate p-values
    p_value_two_tailed = mean(abs(permuted_effects) >= abs(observed_effect));
    p_value_greater = mean(permuted_effects >= observed_effect);
    p_value_less = mean(permuted_effects <= observed_effect);
    
    % Effect size (Cohen's d)
    pooled_std = std([tp1_values; tp2_values]);
    cohens_d = observed_effect / pooled_std;
    
    % Confidence interval (percentile method)
    ci_lower = prctile(permuted_effects, 2.5);
    ci_upper = prctile(permuted_effects, 97.5);
    
    % Package results
    results = struct();
    results.observed_effect = observed_effect;
    results.p_value_two_tailed = p_value_two_tailed;
    results.p_value_greater = p_value_greater;
    results.p_value_less = p_value_less;
    results.cohens_d = cohens_d;
    results.confidence_interval_95 = [ci_lower, ci_upper];
    results.null_distribution = permuted_effects;
    results.n_permutations = nPermutations;
    results.sample_sizes = struct('TP1', length(tp1_values), 'TP2', length(tp2_values));
    results.n_subjects = n_subjects;
    results.subjects = subjects;
    
    % Display results
    % displayResults(results);
    
    % Create plots if requested
    if makePlots
        plotPermutationResults(results);
    end
end

function validateInputData(data)
    % Validate that input data has required fields
    required_fields = {'Subject', 'TimePoint', 'Value'};
    for i = 1:length(required_fields)
        if ~isfield(data, required_fields{i})
            error('Data must contain field: %s', required_fields{i});
        end
    end
    
    % Check that TimePoint only contains 1 and 2
    unique_tp = unique(data.TimePoint);
    if ~all(ismember(unique_tp, [1, 2]))
        error('TimePoint must only contain values 1 and 2');
    end
    
    % Check that all subjects have data for both timepoints
    subjects = unique(data.Subject);
    for s = 1:length(subjects)
        subject_tp = data.TimePoint(data.Subject == subjects(s));
        if ~all(ismember([1, 2], unique(subject_tp)))
            warning('Subject %d does not have data for both timepoints', subjects(s));
        end
    end
end

function displayResults(results)
    % Display formatted results
    fprintf('\n=== PERMUTATION TEST RESULTS ===\n');
    fprintf('Observed Effect (TP2 - TP1): %.4f\n', results.observed_effect);
    fprintf('Cohen''s d: %.4f\n', results.cohens_d);
    fprintf('Two-tailed p-value: %.4f\n', results.p_value_two_tailed);
    fprintf('One-tailed p-value (TP2 > TP1): %.4f\n', results.p_value_greater);
    fprintf('One-tailed p-value (TP2 < TP1): %.4f\n', results.p_value_less);
    fprintf('95%% CI for null distribution: [%.4f, %.4f]\n', ...
            results.confidence_interval_95(1), results.confidence_interval_95(2));
    fprintf('Number of subjects: %d\n', results.n_subjects);
    fprintf('Sample sizes: TP1=%d, TP2=%d\n', ...
            results.sample_sizes.TP1, results.sample_sizes.TP2);
    fprintf('Number of permutations: %d\n', results.n_permutations);
end

function plotPermutationResults(results)
    % Create visualization of permutation test results
    
    figure('Position', [100, 100, 1200, 500]);
    
    % Subplot 1: Histogram of null distribution
    subplot(1, 2, 1);
    histogram(results.null_distribution, 50, 'Normalization', 'probability', ...
              'FaceColor', [0.7, 0.7, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    
    % Add observed effect line
    y_limits = ylim;
    plot([results.observed_effect, results.observed_effect], y_limits, 'r--', ...
         'LineWidth', 2, 'DisplayName', sprintf('Observed Effect: %.4f', results.observed_effect));
    
    % Add zero line
    plot([0, 0], y_limits, 'k-', 'DisplayName', 'Null Hypothesis');
    
    xlabel('Effect Size (TP2 - TP1)');
    ylabel('Probability');
    title('Null Distribution vs Observed Effect');
    legend('Location', 'best');
    grid on;
    grid minor;
    
    % Subplot 2: Q-Q plot
    subplot(1, 2, 2);
    qqplot(results.null_distribution);
    title('Q-Q Plot: Null Distribution vs Normal');
    grid on;
    grid minor;
    
    % Add overall title
    sgtitle(sprintf('Permutation Test Results (p = %.4f)', results.p_value_two_tailed));
end


