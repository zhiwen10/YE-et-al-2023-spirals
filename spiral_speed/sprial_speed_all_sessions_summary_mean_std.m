%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
Fs = 35; % frame sampling rate
scale1 = 8;
padding = 30;
halfpadding = padding/2;
grid_x = [250,350];
grid_y = [350,550];
%%
color_stairs = cbrewer2('seq','BuPu',14);
angular_velocity_all = [];
linear_velocity_all = [];
angular_all = [];
distance_all = [];
angular_velocity_all_all = [];
linear_velocity_all_all = []; 
count1 = 1;
for radiusi = 40:10:100
    for kk = 1:15
    % for kk = 1
        clear filteredSpirals2 filteredSpirals3 tracePhase_raw tracePhase_padded...
            angle_offset_all distance_offset_all  
        %%
        mn = T.MouseID{session_all(kk)};
        tda = T.date(session_all(kk));
        en = T.folder(session_all(kk));
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        serverRoot = expPath(mn, td, en);
        %%
        folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_speed\all_sessions_fftn';
        fname = [mn '_' tdb '_' num2str(en) '.mat'];
        load(fullfile(folder,fname));
        %%
        folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
        fname = [mn '_' tdb '_' num2str(en)];
        dfolder = fullfile(folder,mn,td,num2str(en));
        spiralGroupName = dir(fullfile(dfolder,'*spirals_group_fftn.mat')).name;
        load(fullfile(dfolder,spiralGroupName));
        clear index1 indx2 groupedCells
        indx2 = cellfun(@(x) size(x,1), archiveCell);
        groupedCells = archiveCell(indx2>=2);
        filteredSpirals2 = cell2mat(groupedCells);
        %%
        index1 = (filteredSpirals2(:,2)>grid_x(1) & filteredSpirals2(:,2)<grid_x(2)...
            &filteredSpirals2(:,1)>grid_y(1) & filteredSpirals2(:,1)<grid_y(2));
        filteredSpirals3 = filteredSpirals2(index1,:);
        %%
        % indx2 = cellfun(@(x) size(x,2), distance_offset_all);
        indx2 = (filteredSpirals3(:,3) == radiusi);
        %%
        angle_offset_all1 = cellfun(@(x) mean(x,'omitnan'),angle_offset_all,'UniformOutput',false);
        angle_offset_all_select = angle_offset_all1(indx2);
        angle_offset_all_sessions{kk} = angle_offset_all_select;    
        %%
        distance_offset_all_select = distance_offset_all(indx2); % radius of 8*10=80 pixels
        distance_offset_all_sessions{kk} = distance_offset_all_select;
    end
    %
    angle_offset_all_sessions2 = horzcat(angle_offset_all_sessions{:});
    distance_offset_all_sessions2 = horzcat(distance_offset_all_sessions{:});
    angular_velocity_all = vertcat(angle_offset_all_sessions2{:})*35;
    distance_offset_all1 = cellfun(@(x) mean(x,1,'omitnan'),distance_offset_all_sessions2,'UniformOutput',false);
    linear_velocity_all= vertcat(distance_offset_all1{:});
%     mean_angle = mean(angular_velocity_all);
%     std_angle = std(angular_velocity_all);
%     mean_linear = mean(linear_velocity_all);
%     std_linear = std(linear_velocity_all);
    %%
    angular_velocity_all_max = angular_velocity_all(:,end);
    max_angle = mean(angular_velocity_all_max);
    std_angle = std(angular_velocity_all_max);
    sem_angle = std_angle/sqrt(numel(angular_velocity_all_max));
    angular_all(:,count1) = [max_angle,std_angle,sem_angle];
    
    linear_velocity_all_max = linear_velocity_all(:,end);
    max_distance = mean(linear_velocity_all_max);
    std_distance = std(linear_velocity_all_max);
    sem_distance = std_distance/sqrt(numel(linear_velocity_all_max));
    distance_all(:,count1) = [max_distance,std_distance, sem_distance];
    %%
    angular_velocity_all_all = [angular_velocity_all_all;angular_velocity_all_max];
    linear_velocity_all_all = [linear_velocity_all_all;linear_velocity_all_max];
    
    count1 = count1+1;
end
%%
mean_anglular_all = mean(angular_velocity_all_all);
std_angular_all = std(angular_velocity_all_all);
%%
frq = angular_velocity_all_all/(2*pi);
mean_frq = mean(frq);
std_frq = std(frq);
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 500 300])
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
print(h1,'spiral_speed_all_mean+std', '-dpdf', '-bestfit', '-painters');
save('spiral_speed_radius_mean_std.mat','angular_all','distance_all');