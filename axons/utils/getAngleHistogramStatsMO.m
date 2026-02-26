function [angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStatsMO(T1)
scale = 1;
scale1 = 15;
center = [377,428]; % MOp
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
orthog_vector = soma_center1-center;
angle_diff_real = deg2rad(T1.center_bias_angle)';
angle_diff_perm = [];
for k = 1:100
    rp_indx = randperm(size(T1.axon_bias_1,1));
    vector_all1 = [T1.axon_bias_1(rp_indx,:),T1.axon_bias_2(rp_indx,:)];
    for i =1:size(T1.axon_bias_1,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
    angle_diff_perm(k,:) = ang_diff1;
end