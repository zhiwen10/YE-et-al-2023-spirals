function getSpectrumSpiralRegression(data_folder,save_folder)
load(fullfile(data_folder, 'spirals\spectrum','spirals_fft_radius.mat'));
powerV1 = readNPY(fullfile(data_folder,'spirals\spectrum','powerSpectrum.npy'));
freq = readNPY(fullfile(data_folder,'spirals\spectrum','frequency.npy'));
color1 = cbrewer2('seq','Reds',15);
%%
spirals_sum = sum(count_sample_control(:,4:10),2);
%%
area_names = {'VISp','RSP','SSp-ul','SSp-ll','SSp-m','SSp-n','SSp-bfd'};
%%
exponents = zeros(7,15);
center_freq = zeros(7,15);
power = zeros(7,15);
bandwidth = zeros(7,15);
offset = zeros(7,15);
for i = 1:15
    data_folder1 = fullfile(data_folder,'spirals\spectrum\fooof');
    Ta = readtable(fullfile(data_folder1,['file' num2str(i-1) '.csv']));
    exponents(:,i) = Ta.exponent;
    center_freq(:,i) = Ta.cf_0;
    power(:,i) = Ta.pw_0;
    bandwidth(:,i) = Ta.bw_0;
    offset(:,i) = Ta.offset;
end
power(isnan(power)) = 0;
%%
offset1z = zscore(offset(6,:),[],2);
exponents1z = zscore(exponents(6,:),[],2);
power1z = zscore(power(6,:),[],2);
%%
regressor = [offset1z;exponents1z;power1z]';
regressor1 = [regressor,ones(size(regressor,1),1)];
%%
weights = regressor1\spirals_sum;
%% unique variance explained
mdl = fitlm(regressor,spirals_sum);
r2_total = mdl.Rsquared.Ordinary;
% leave one out regression
regressor_index = [2,3;1,3;1,2];
for i = 1:3
    regressor2 = regressor(:,regressor_index(i,:));
    mdl = fitlm(regressor2,spirals_sum);
    r2(i,1) = mdl.Rsquared.Ordinary;
end
r2_unique = r2_total-r2;
%% exponents variance
mdl = fitlm(exponents1z',spirals_sum);
r2_exponents = mdl.Rsquared.Ordinary;