function plotSpiralvsFooof(data_folder,save_folder)
load(fullfile(data_folder, 'spirals\spectrum','spirals_fft_radius.mat'));
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
exponents = nan(7,15);
center_freq = nan(7,15);
power = nan(7,15);
bandwidth = nan(7,15);
for i = 1:15
    data_folder1 = fullfile(data_folder,'spirals\spectrum\fooof');
    Ta = readtable(fullfile(data_folder1,['file' num2str(i-1) '.csv']));
    exponents(:,i) = Ta.exponent;
    center_freq(:,i) = Ta.cf_0;
    power(:,i) = Ta.pw_0;
    bandwidth(:,i) = Ta.bw_0;
end
%%
mean_percentage_change = mean(percentage_change,1);
std_percentage_change = std(percentage_change,[],1);
sem_percentage_change = std_percentage_change/sqrt(15);
%%
x_point = normrnd(0,0.1,[1,15]);
scatter_x = repmat([1:10],15,1);
scatter_x = scatter_x+repmat(x_point',1,10);
%%
% [A,I] = sort(count_sample_control(:,10));
% count_sample_control = count_sample_control(I,:);
% count_sample_permute = count_sample_permute(I,:);
% percentage_change = percentage_change(I,:);
%%
color1 = cbrewer2('seq','YlOrRd',15);
h1 = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
for radius = 1:10
    subplot(1,21,radius)
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    scatter(1,control_all,4,color1,'filled');
    hold on;
    scatter(2,permute_all,4,color1,'filled');
    hold on;
    % x_control = ones(size(control_all)); x_permute = ones(size(control_all))*2;
    for kk = 1:15
        plot([1,2]',[control_all(kk),permute_all(kk)]','Color',color1(kk,:));
    end
    ylim([0,2.5]);
end
subplot(1,21,[12:21]);
for radius = 1:10
    scatter(scatter_x(:,radius),percentage_change(:,radius)*100,6,color1,'filled');
    hold on;
end
x = 1:10;
bar(x,mean_percentage_change*100,'FaceColor','None');
hold on
er = errorbar(x,mean_percentage_change*100,sem_percentage_change*100);    
er.Color = [0,0,0];                         
er.LineStyle = 'none';  
hold off
ylim([-100 500]);
xlabel('Radius');
ylabel('Percentage change (%)');
%%
% bandwidth = bandwidth(:,I);
% center_freq = center_freq(:,I);
% power = power(:,I);
% exponents = exponents(:,I);
power(isnan(power))=0;
%%
figure;
for radius = 1:6
    for i = 1:7
        subplot(6,7,i+(radius-1)*7);
        % scatter(count_sample_control(:,radius+4),power(i,:),6,color1,'filled');
        scatter(count_sample_control(:,radius+4),power(i,:),6,'k','filled');
    %     xlim([1,3]);
    %     ylim([0,1]);
        title(area_names{i});
    end
end
%%
figure;
for i = 1:6
    subplot(1,6,i);
    scatter(power(1,:),power(6,:),12,count_sample_control(:,i+4),'filled');
    xlabel('VISp power');
end
ylabel('SSpn power');
colorbar;
%%

