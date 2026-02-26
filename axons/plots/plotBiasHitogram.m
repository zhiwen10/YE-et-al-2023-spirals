function h2d = plotBiasHitogram(data_folder,save_folder)
%%
T = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells.csv'));
h2d = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
ax1 =subplot(1,1,1);
soma_center1 = [T.soma_center_1,T.soma_center_2];
axon_bias = [T.axon_bias_1,T.axon_bias_2];
edges = 0:10:180;
[N1,edges] = histcounts(T.center_bias_angle,edges);   
mean_bias = mean(T.center_bias_angle);
% shuffle
center = [244,542]; % SSp-un
% center = [236,629]; % SSp-un
orthog_vector = soma_center1-center;
for k = 1:100
    rp_indx = randperm(size(axon_bias,1));
    vector_all1 = axon_bias(rp_indx,:);
    for i =1:size(axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
ax1 = plotAngleStairsShuffle(ax1,edges,N1,N_rp);
hold on;
xline(mean_bias,'r');
xlim([0,180]); 
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};
%%
[angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStats(T);
[h mu ul ll] = circ_mtest(angle_diff_real,deg2rad(90));
mean_rad = [mu,ll,ul];
mean_deg = rad2deg(mean_rad);
% for i = 1:100
for i = 1
    [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end
%%
print(h2d, fullfile(save_folder,'Fig2d_axon_bias_histogram.pdf'),...
    '-dpdf', '-bestfit', '-painters');