function [hs3f,hs3g] = plotScatterDataVsFftn2(T,data_folder,save_folder,freq)
%%
freq_name = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
local_data_folder = fullfile(data_folder,...
    'spirals','spirals_fftn');
control_data_folder = fullfile(local_data_folder,...
    freq_name,'control_stats');
fftn_data_folder = fullfile(local_data_folder,...
    freq_name,'fftn_stats');
%%
color1 = 'k';
ssp_index = get_ssp_index(data_folder);
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
count_sample_control = zeros(15,10);
count_sample_permute = zeros(15,10);
for kk = 1:size(T,1)
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    control = load(fullfile(control_data_folder,[fname '_density.mat']));
    permute = load(fullfile(fftn_data_folder,[fname '_density.mat']));
    count = 1;
    for radius = 10:10:100
        clear unique_spirals_unit_control unique_spirals_unit_permute 
        spirals_control = control.spiral_density{count};
        spirals_permute = permute.spiral_density{count};
        [lia,locb] = ismember(spirals_control(:,1:2),ssp_index,'rows');
        spirals_control = spirals_control(lia,:);
        [lia1,locb1] = ismember(spirals_permute(:,1:2),ssp_index,'rows');
        spirals_permute = spirals_permute(lia1,:);
        unique_spirals_unit_control = spirals_control(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = spirals_permute(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        % interp histgram counts 
        count_sample_control(kk,count) = max(unique_spirals_unit_control(:));
        count_sample_permute(kk,count) = max(unique_spirals_unit_permute(:));
        count = count+1;
    end
end
%%
hs3f = figure;
for radius = 1:10
    subplot(1,10,radius)
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    
    x_control = ones(size(control_all)); 
    x_permute = ones(size(control_all))*2;
    scatter(x_control,control_all,4,color1,'filled');
    hold on;
    scatter(x_permute,permute_all,4,color1,'filled');
    hold on;

%     plot([x_control,x_permute]',[control_all,permute_all]','k');
    for kk = 1:15
        % plot([1,2]',[control_all(kk),permute_all(kk)]','Color',color1(kk,:));
        plot([1,2]',[control_all(kk),permute_all(kk)]','Color',color1);
    end
    ylim([0,2.5]);
end
%%
print(hs3f, fullfile(save_folder,...
    ['FigS3f_spirals_across_session_' freq_name '.pdf']),...
    '-dpdf', '-bestfit', '-painters');
%%
for radius = 1:10
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    [ha(radius),pa(radius)] = ttest(control_all,permute_all);
end
%% percentage change
for radius = 1:10
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    percentage_change(:,radius) = (control_all-permute_all)./permute_all;
end
%%
for i = 1:10
    percentage_temp = percentage_change(:,i);
    [h3(i),p3(i)] = ttest(percentage_temp);
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
hs3g = figure;
% scatter(scatter_x,percentage_change*100,6,'k','filled');
for radius = 1:10
    scatter(scatter_x(:,radius),percentage_change(:,radius)*100,6,color1,'filled');
    hold on;
end
hold on;
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
print(hs3g, fullfile(save_folder,...
    ['FigS3g_spirals_change_across_session_' freq_name '.pdf']), ...
    '-dpdf', '-bestfit', '-painters');
end