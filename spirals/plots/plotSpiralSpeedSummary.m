%%
function h1fg = plotSpiralSpeedSummary(T,data_folder,save_folder)
load(fullfile(data_folder,'spirals','spirals_speed','speed_all.mat'));
%%
angular_all = [];
distance_all = [];
count1 = 1;
for radiusi = 40:10:100
    clear angle_offset_all distance_offset_all
    angle_offset_all = horzcat(angle_offset{:,count1});
    distance_offset_all = horzcat(distance_offset{:,count1});
    %%
    angular_velocity_temp = vertcat(angle_offset_all{:})*35;
    distance_offset_all1 = cellfun(@(x) mean(x,1,'omitnan'),...
        distance_offset_all,'UniformOutput',false);
    linear_velocity_temp= vertcat(distance_offset_all1{:});
    %%
    angular_velocity_all_max = angular_velocity_temp(:,end);
    max_angle = mean(angular_velocity_all_max);
    std_angle = std(angular_velocity_all_max);
    sem_angle = std_angle/sqrt(numel(angular_velocity_all_max));
    angular_all(:,count1) = [max_angle,std_angle,sem_angle];

    linear_velocity_all_max = linear_velocity_temp(:,end);
    max_distance = mean(linear_velocity_all_max);
    std_distance = std(linear_velocity_all_max);
    sem_distance = std_distance/sqrt(numel(linear_velocity_all_max));
    distance_all(:,count1) = [max_distance,std_distance, sem_distance];
    %%
    count1 = count1+1;
end
%%
h1fg = figure('Renderer', 'painters', 'Position', [100 100 500 300])
pixSize = 3.45/1000/0.6*3;  % mm / pix
x_value = [40:10:100]*pixSize;
subplot(1,2,1);
errorbar(x_value,angular_all(1,:), angular_all(2,:),'r',"CapSize",8);
xlim([0.5,2]); ylim([0,80]);
t1 = [0.5, 1.0, 1.5, 2.0];
t2 = cellstr(string(t1));
xticks(t1); xticklabels(t2);
xlabel('Spiral radius (mm)');
ylabel('Angular speed (rad/s)');
subplot(1,2,2);
errorbar(x_value,distance_all(1,:),distance_all(2,:),'r',"CapSize",8);
xlim([0.5,2]); ylim([0,80]);
t1 = [0.5, 1.0, 1.5, 2.0];
t2 = cellstr(string(t1));
xticks(t1); xticklabels(t2);
xlabel('Spiral radius (mm)');
ylabel('Linear speed (mm/s)');
%%
print(h1fg,fullfile(save_folder,'Fig1fg_spiral_speed_all.pdf'),...
    '-dpdf', '-bestfit', '-painters');
