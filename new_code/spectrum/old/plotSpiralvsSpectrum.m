function plotSpiralvsSpectrum(data_folder,save_folder)
load(fullfile(data_folder, 'spirals\spectrum','spirals_fft_radius.mat'));
powerV1 = readNPY(fullfile(data_folder,'spirals\spectrum','powerSpectrum.npy'));
freq = readNPY(fullfile(data_folder,'spirals\spectrum','frequency.npy'));
color1 = cbrewer2('seq','Reds',15);
%%
spirals_sum = sum(count_sample_control(:,4:10),2);
%%
[A,I] = sort(spirals_sum);
% [A,I] = sort(count_sample_control(:,8));
powerV1 = powerV1(:,:,I);
%%
pixels(1,:) = [845,835]; % VISp
pixels(2,:) = [775,650]; % RSP
pixels(3,:) = [590,750]; % SSp-ul
pixels(4,:) = [520,850]; % SSp-ll
pixels(5,:) = [480,950]; % SSp-m
pixels(6,:) = [550,950]; % SSp-n
pixels(7,:) = [675,905]; % SSp-bfd
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
offset = offset(:,I);
exponents = exponents(:,I); 
power = power(:,I);
power(isnan(power))=0;
center_freq = center_freq(:,I);
%%
color2 = cbrewer2('seq','YlOrRd',9);
scale3 = 5; hemi = []; lineColor = 'w';
h1 = figure('Renderer', 'painters', 'Position', [100 100 1200 700]);
xticks_value = log10([0.5,2,4,6,8]);
for kk = 1:15
    for i = 1:7
        ax(i) = subplot(5,7,i);
        powerV = log10(squeeze(powerV1(:,i,kk))); 
        plot(log10(freq(2:end)),powerV(2:end),'color',color1(kk,:));
        xlabel('Frequency (Hz)'); ylabel('log10(Power)');
        xticks(xticks_value);
        xticklabels({'0.5','2','4','6','8'});
        xlim([xticks_value(1),xticks_value(5)]);
        title(area_names{i});
        hold on;  
    end
end
for i = 1:7
    %%
    ax1(i) = subplot(5,7,i+7);
    scatter(A,offset(i,:),6,color1,'filled');
    ylim([-5.5,-4]);
    xlabel('Spirals density'); ylabel('Offset');
    %%
    ax2(i) = subplot(5,7,i+14);
    scatter(A,exponents(i,:),6,color1,'filled');
    ylim([1,3]);
    xlabel('Spirals density'); ylabel('Exponents');
    %%
    ax3(i) = subplot(5,7,i+21);
    scatter(A,center_freq(i,:),6,color1,'filled');
    ylim([0,6]);
    xlabel('Spirals density'); ylabel('CenterFreq');
    %%
    ax4(i) = subplot(5,7,i+28);
    scatter(A,power(i,:),6,color1,'filled');
    ylim([0,0.6]);
    xlabel('Spirals density'); ylabel('Power');
end
%%
print(h1, fullfile(save_folder,'power_spectrum_params'), '-dpdf', '-bestfit', '-painters');