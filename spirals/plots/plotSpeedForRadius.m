function h9gj = plotSpeedForRadius(radius,data_folder,save_folder)
%%
load(fullfile(data_folder,'spirals','spirals_speed',...
    ['speed_' num2str(radius) 'pixels.mat']));
%%
mean_angle = mean(angular_velocity_all);
std_angle = std(angular_velocity_all);
mean_linear = mean(linear_velocity_all);
std_linear = std(linear_velocity_all);
pixSize = 3.45/1000/0.6*3;  % mm / pix
radius2 = round(radius*pixSize*10);
%%
spiral_n = size(angular_velocity_all,1);
angular_velocity_per_spiral = nan(spiral_n,12); % maximum length is 12
%%
h9gj = figure('Renderer', 'painters', 'Position', [100 100 1100 250]);
ax1 = subplot(1,4,1);
p = randperm(spiral_n,1000);
scale1 = 8;
pixSize = 3.45/1000/0.6*3;  % mm / pix
sampling_radius = [1:size(angular_velocity_all,2)]*scale1*pixSize;
sampling_radius_all = [1:12]*scale1*pixSize;
for i = 1:1000
    current_linear_velocity1(i,:) = angular_velocity_all(p(i),:);
    sampling_radius1(i,:) = sampling_radius;
end
scatter_kde(sampling_radius1(:), current_linear_velocity1(:), 'filled', 'MarkerSize', 6);
colormap(ax1,parula);

mean_angle_offset_all = mean(angular_velocity_all,1,'omitnan');
hold on;
plot(sampling_radius,mean_angle_offset_all,'-o','color',[1, 0,0],...
        'MarkerSize',3,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none');
xlim([0,1.6]);
ylim([0,100]);
yticks([0 20 40 60 80 100]);
yticklabels({'0' '20' '40' '60' '80' '100'});
xticks([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6]);
xticklabels({'0' '0.2' '0.4' '0.6' '0.8' '1.0' '1.2' '1.4' '1.6'});
xlabel('Radius (mm)');
ylabel('Angular velocity (rad/s)');

subplot(1,4,2);
edges = 0:2:100;
angular_velocity_all_max = angular_velocity_all(:,end);
[N1,edges1] = histcounts(angular_velocity_all_max,edges);
max_angle = mean(angular_velocity_all_max);
std_angle = std(angular_velocity_all_max);
a1 = stairs(edges1(1:end-1), N1);
hold on;
xline(max_angle,'r');
text(max_angle,10,num2str(round(max_angle*10)/10));
xlim([0,100]);
xlabel('Angular velocity (rad/s)');
ylabel('spiral counts');

ax3 = subplot(1,4,3);
distance_offset_per_spiral = size(linear_velocity_all,1);
radius_n = size(linear_velocity_all,2);
for i = 1:1000
    current_linear_velocity1(i,:) = linear_velocity_all(p(i),:);
    sampling_radius1(i,:) = sampling_radius;
end
scatter_kde(sampling_radius1(:), current_linear_velocity1(:), 'filled', 'MarkerSize', 6);
colormap(ax3,parula);

mean_distance_offset_all = mean(linear_velocity_all,1,'omitnan');
hold on;
plot(sampling_radius,mean_distance_offset_all,'-o','color',[1, 0,0],...
        'MarkerSize',3,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none');
xlim([0,1.6]);
ylim([0,100]);
yticks([0 20 40 60 80 100]);
yticklabels({'0' '20' '40' '60' '80' '100'});
xticks([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6]);
xticklabels({'0' '0.2' '0.4' '0.6' '0.8' '1.0' '1.2' '1.4' '1.6'});
xlabel('Radius (mm)');
ylabel('Speed(mm/s)');

%
subplot(1,4,4)
linear_velocity_all_max = linear_velocity_all(:,end);
max_distance = mean(linear_velocity_all_max);
std_distance = std(linear_velocity_all_max);
edges2 = edges*radius*pixSize;                                             % convert from angular velocity edges to linear
a2 = stairs(edges2(1:end-1), N1);
hold on;
xline(max_distance,'r');
xlim([0,100]);
text(max_distance,10,num2str(round(max_distance*10)/10));
xlabel('Speed(mm/s)');
ylabel('spiral counts');
%%
print(h9gj,fullfile(save_folder,...
    ['Figs9gj_spiralSpeed_' num2str(radius) 'pixels']),...
    '-dpdf', '-bestfit', '-painters');
