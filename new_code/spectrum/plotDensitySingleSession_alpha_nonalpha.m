%%
save_folder1 = fullfile(data_folder,'spirals\spirals_freq\spirals_fft_bump_density\singleSessions');
spirals_folder = fullfile(data_folder,'spirals\spirals_freq\spirals_fft_bump');
label1 = 'alpha'; label2 = 'nonalpha';
frame_all = 0;
hist_bin = 40;
for kk = 1:size(T,1)
    %% session info
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    % load data
    fname = [mn '_' tdb '_' num2str(en)];     
    load(fullfile(spirals_folder,[fname '_histogram.mat']));
    %%
    clear spiral_density
    count = 1;
    for radius = 40:10:100
        clear unique_spirals
        spirals_temp = [];
        spirals_temp = filteredSpirals_alpha(filteredSpirals_alpha(:,3)==radius,:);
        [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
            density_color_plot(spirals_temp,hist_bin);
        spiral_density{count} = unique_spirals;
        count = count+1;
    end
    frameN = numel(alphaFrames);
    save(fullfile(save_folder1,label1,[fname '_density.mat']), ...
        'spiral_density','frameN');
    %%
    clear spiral_density
    count = 1;
    for radius = 40:10:100
        clear unique_spirals
        spirals_temp = [];
        spirals_temp = filteredSpirals_nonalpha(filteredSpirals_nonalpha(:,3)==radius,:);
        [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
            density_color_plot(spirals_temp,hist_bin);
        spiral_density{count} = unique_spirals;
        count = count+1;
    end
    frameN = numel(nonalphaFrames);
    save(fullfile(save_folder1,label2,[fname '_density.mat']), ...
        'spiral_density','frameN');
end
%%
local_data_folder = fullfile(data_folder,...
    'spirals\spirals_freq\spirals_fft_bump_density\singleSessions');
alpha_data_folder = fullfile(local_data_folder,'alpha');
nonalpha_data_folder = fullfile(local_data_folder,'nonalpha');
color1 = 'k';
ssp_index = get_ssp_index(data_folder);
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
count_sample_alpha = zeros(15,7);
count_sample_nonalpha = zeros(15,7);
for kk = 1:size(T,1)
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    alpha = load(fullfile(alpha_data_folder,[fname '_density.mat']));
    nonalpha = load(fullfile(nonalpha_data_folder,[fname '_density.mat']));
    count = 1;
    for radius = 40:10:100
        clear unique_spirals_unit_control unique_spirals_unit_permute 
        spirals_alpha = alpha.spiral_density{count};
        spirals_nonalpha = nonalpha.spiral_density{count};
        [lia,locb] = ismember(spirals_alpha(:,1:2),ssp_index,'rows');
        spirals_alpha = spirals_alpha(lia,:);
        [lia1,locb1] = ismember(spirals_nonalpha(:,1:2),ssp_index,'rows');
        spirals_nonalpha = spirals_nonalpha(lia1,:);
        unique_spirals_unit_alpha = spirals_alpha(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_alpha = unique_spirals_unit_alpha./alpha.frameN*35; % spirals/(mm^2*s)
        unique_spirals_unit_nonalpha = spirals_nonalpha(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_nonalpha = unique_spirals_unit_nonalpha./nonalpha.frameN*35; % spirals/(mm^2*s)
        % interp histgram counts 
        count_sample_alpha(kk,count) = max(unique_spirals_unit_alpha(:));
        count_sample_nonalpha(kk,count) = max(unique_spirals_unit_nonalpha(:));
        count = count+1;
    end
end
%%
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_fft_bump_density');
hs3d = figure;
for radius = 1:7
    subplot(1,7,radius)
    alpha_all = squeeze(count_sample_alpha(:,radius));
    nonalpha_all = squeeze(count_sample_nonalpha(:,radius));
    scatter(1,alpha_all,4,color1,'filled');
    hold on;
    scatter(2,nonalpha_all,4,color1,'filled');
    hold on;
%     x_control = ones(size(control_all)); x_permute = ones(size(control_all))*2;
%     plot([x_control,x_permute]',[control_all,permute_all]','k');
    for kk = 1:15
        % plot([1,2]',[control_all(kk),permute_all(kk)]','Color',color1(kk,:));
        plot([1,2]',[alpha_all(kk),nonalpha_all(kk)]','Color',color1);
    end
    ylim([0,2.5]);
end
print(hs3d, fullfile(save_folder,...
    ['Spirals_across_session.pdf']),'-dpdf', '-bestfit', '-painters');
%%
for radius = 1:7
    alpha_all = squeeze(count_sample_alpha(:,radius));
    nonalpha_all = squeeze(count_sample_nonalpha(:,radius));
    [ha(radius),pa(radius)] = ttest(alpha_all,nonalpha_all);
end
%% percentage change
for radius = 1:7
    alpha_all = squeeze(count_sample_alpha(:,radius));
    nonalpha_all = squeeze(count_sample_nonalpha(:,radius));
    percentage_change(:,radius) = (alpha_all-nonalpha_all)./nonalpha_all;
end
%%
for i = 1:7
    percentage_temp = percentage_change(:,i);
    [h4(i),p4(i)] = ttest(percentage_temp);
end
%%
mean_percentage_change = mean(percentage_change,1);
std_percentage_change = std(percentage_change,[],1);
sem_percentage_change = std_percentage_change/sqrt(15);
%%
x_point = normrnd(0,0.1,[1,15]);
scatter_x = repmat([1:7],15,1);
scatter_x = scatter_x+repmat(x_point',1,7);
%%
hs3e = figure;
% scatter(scatter_x,percentage_change*100,6,'k','filled');
for radius = 1:7
    scatter(scatter_x(:,radius),percentage_change(:,radius)*100,6,color1,'filled');
    hold on;
end
hold on;
x = 1:7;
bar(x,mean_percentage_change*100,'FaceColor','None');
hold on
er = errorbar(x,mean_percentage_change*100,sem_percentage_change*100);    
er.Color = [0,0,0];                         
er.LineStyle = 'none';  
hold off
ylim([-100 500]);
xlabel('Radius');
ylabel('Percentage change (%)');
print(hs3e, fullfile(save_folder,...
    'Spirals_change_across_session.pdf'),'-dpdf', '-bestfit', '-painters');
