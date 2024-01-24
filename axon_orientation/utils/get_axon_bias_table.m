function T = get_axon_bias_table(soma_center,vector_all,angle1,polarity,labels)
%% calculate angular position of soma
clear ang_diff1 ang_diff2
center = [244,542]; % SSp-un
orthog_vector = soma_center-center;
for i =1:size(vector_all,1)
    u = [vector_all(i,:) 0];
    v = [orthog_vector(i,:) 0];
    ang_diff1(i,1) = atan2(norm(cross(u,v)),dot(u,v));
end    
ang_diff2 = rad2deg(ang_diff1);
%% calculate angular position of axon
r = 1;
soma_angle = atan2(orthog_vector(:,2),orthog_vector(:,1));
z = r*exp(1i*soma_angle);
%% sort neuron by soma position
T = table(soma_center,soma_angle,z,vector_all,angle1,polarity,ang_diff2,labels);
T = renamevars(T,["z","vector_all","angle1","polarity","ang_diff2"], ...
                 ["soma_polar_angle","axon_bias","bias_angle","pc_ratio","center_bias_angle"]);
% writetable(T,"axon_bias_table.csv");
