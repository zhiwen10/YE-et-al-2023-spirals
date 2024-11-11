%%
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO2.csv'));
%%
axon_vector = [T1.axon_bias_1,T1.axon_bias_2];
theta = 90;
for i = 1:size(axon_vector,1)
    v1 = axon_vector(i,:);
    vr = v1*[1;1i]*exp(-1i*theta*pi/180);
    angle1(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
end
angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
%%
center = [377,428]; % MOp
angle2 = angle1(T1.soma_center_1<center(1));
angle3= angle1(T1.soma_center_1>center(1));
%%
load(fullfile(data_folder,'Revision','axons','SSp_roi.mat'));
T3 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
tf = inROI(roi,T3.soma_center_1,T3.soma_center_2);
T3 = T3(tf,:);
axon_vector = [T3.axon_bias_1,T3.axon_bias_2];
theta = 90;
for i = 1:size(axon_vector,1)
    v1 = axon_vector(i,:);
    vr = v1*[1;1i]*exp(-1i*theta*pi/180);
    angle4(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
end
angle4(angle4>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
%%
figure;
subplot(1,3,1);
histogram(angle2,20);
xlim([-90,90]);
ylim([0,20]);
subplot(1,3,2);
histogram(angle3,20);
xlim([-90,90]);
ylim([0,20]);
subplot(1,3,3);
histogram(angle4,20);
xlim([-90,90]);
ylim([0,20]);