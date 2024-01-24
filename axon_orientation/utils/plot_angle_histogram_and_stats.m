function [ax4,angle_diff_real,angle_diff_perm,edges,N1,N_rp] = plot_angle_histogram_and_stats(ax4,T1,coords)
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
center = [244,542]; % SSp-un
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T1.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = T1.soma_center;
cell_n = size(T1,1);
scale1 = 15;
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
center = [244,542]; % SSp-un
orthog_vector = soma_center1-center;
angle_diff_real = deg2rad(T1.center_bias_angle)';
angle_diff_perm = [];
for k = 1:100
    rp_indx = randperm(size(T1.axon_bias,1));
    vector_all1 = T1.axon_bias(rp_indx,:);
    for i =1:size(T1.axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
    angle_diff_perm(k,:) = ang_diff1;
end
ax4 = plot_angle_stairs_shuffle(ax4,edges,N1,N_rp);
hold on;
xline(mean_bias,'r');
xlim([0,180]); 
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};
% ylim([0,40]);