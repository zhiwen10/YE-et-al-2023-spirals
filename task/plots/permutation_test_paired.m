%% compare correct and incorrect for post response
T = readtable('spirals_pre_post2.csv');
T = T(T.TimePoint == 2,:); % only look at post response
T = T(:,[1,2,4,5]);
data.Subject = T.Subject; % data structure
data.TimePoint = T.Outcome; % rename outcome to timepoints for function variable convinience
data.Value = T.Value;
data.Sessions = T.Session;
results = permutation_test_paired_timepoints(data);
%% compare correct and incorrect for post-pre 
T = readtable('spirals_pre_post2.csv');
diff = T{T.TimePoint == 2,5}-T{T.TimePoint == 1,5}; % post-pre
T = T(T.TimePoint == 2,:); %use labels at post response
T = T(:,[1,2,4,5]);
data.Subject = T.Subject; % data structure
data.TimePoint = T.Outcome; % rename outcome to timepoints for function variable convinience
data.Value = diff; % use diff as value
data.Sessions = T.Session;
results = permutation_test_paired_timepoints(data);
