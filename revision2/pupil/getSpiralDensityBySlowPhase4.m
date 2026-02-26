function getSpiralDensityBySlowPhase2(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01;                                                            % mm/pix
pixArea = pixSize^2;
ssp_index = get_ssp_index(data_folder);
hist_bin = 40;
%%
load('spirals_sort_by_slow_phase3.mat');
%%
mouseN = size(T,1);
count_sample = zeros(mouseN,size(spirals_sort,2));
for kk  = 1:mouseN
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    t = readNPY(fullfile(session_root,'svdTemporalComponents_corr.timestamps.npy'));
    %%
    for m = 1:size(spirals_sort,2)
        clear spirals_temp unique_spirals unique_spirals_unit
        spirals_temp = spirals_sort{kk,m};
        [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
        [lia,locb] = ismember(unique_spirals(:,1:2),ssp_index,'rows');
        unique_spirals = unique_spirals(lia,:);
        unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit = unique_spirals_unit./numel(t)*35*size(spirals_sort,2);                % spirals/(mm^2*s)
        count_sample(kk,m) = max(unique_spirals_unit(:));
    end
end
%%
pairs = [1,2];
for i = 1:size(pairs,1)
    [p1(i),h1(i)] = ttest(count_sample(:,pairs(i,1)),count_sample(:,pairs(i,2)));
end
%%
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
h1= figure('Renderer', 'painters', 'Position', [100 100 300 500]);
subplot(1,1,1);
for phaseBin = 1:size(spirals_sort,2)
    control_all = squeeze(count_sample(:,phaseBin));    
    x_control = ones(size(control_all))*phase_centers(phaseBin);     
    scatter(x_control,control_all,4,'k','filled');
    hold on;
end
pairs2 = [1:n_bins-1;2:n_bins]';
for i = 1:size(pairs2,1)
    control_all1 = squeeze(count_sample(:,pairs2(i,1)));  
    control_all2 = squeeze(count_sample(:,pairs2(i,2)));  
    for kk = 1:mouseN       
        plot(phase_centers(pairs2(i,:))',[control_all1(kk),control_all2(kk)]','k');
        hold on;
    end
end
xticks([-pi:pi/2:pi]);
xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
xlim([-pi,pi]);
% xticks([1:6]);
% xticklabels({'Up','Down'});
ylim([0,8]);
yticks([0:1:8]);
yticklabels([0:1:8]);
ylabel('Peak spiral density (spirals/mm2*s)');
count_mean = mean(count_sample,1);
count_sem = std(count_sample,[],1)./sqrt(size(count_sample,1));
hold on;
shadedErrorBar(phase_centers, count_mean, count_sem, 'lineprops', '-g')
% hold on;
% errorbar(phase_centers, count_mean, count_sem,'g');
%%
print(h1, 'pupil_phase_spirals.pdf','-dpdf', '-bestfit', '-painters');
%%
ratio_spiral = count_sample(:,1)-count_sample(:,9);
figure;
scatter(ratio,ratio_spiral);