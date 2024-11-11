%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
center = [377,428]; % MOp
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
T3 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
tf = inROI(roi,T3.soma_center_1,T3.soma_center_2);
T3 = T3(tf,:);
axon_vector = [T3.axon_bias_1,T3.axon_bias_2];
theta = 90;
for i = 1:size(axon_vector,1)
    v1 = axon_vector(i,:);
    vr = v1*[1;1i]*exp(-1i*theta*pi/180);
    angle1(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
end
angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
angle2a = angle1(T3.soma_center_1<center(1) & T3.soma_center_2>center(2));
angle2b = angle1(T3.soma_center_1>center(1));
angle2c = angle1(T3.soma_center_1<center(1) & T3.soma_center_2<center(2));
%
figure;
subplot(1,3,1);
hist2a = histogram(angle2a,10,'faceColor','None','EdgeColor','k','DisplayStyle','stairs');
angle2a_max = mean(hist2a.BinEdges(hist2a.Values == max(hist2a.Values)));
xline(angle2a_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
subplot(1,3,2);
hold on;
hist2b = histogram(angle2b,10,'faceColor','None','EdgeColor','r','DisplayStyle','stairs');
angle2b_max = hist2b.BinEdges(hist2b.Values == max(hist2b.Values));
xline(angle2b_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
subplot(1,3,3);
hold on;
hist2c = histogram(angle2c,10,'faceColor','None','EdgeColor','b','DisplayStyle','stairs');
angle2c_max = hist2c.BinEdges(hist2c.Values == max(hist2c.Values));
xline(angle2c_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
%%
angle_all = [angle2a_max+5,angle2b_max+5,angle2c_max+5];
angle_all2 = angle_all+90;
angle_all2(2:3) = angle_all2(2:3)-180;
color2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,angle_all2);
%
angle_all3 = angle_all;
angle_all3(2) = angle_all3(2)-90;
r = 1;
z = r*exp(1i*deg2rad(angle_all3));
z=z';
%
index2a = (T3.soma_center_1<center(1) & T3.soma_center_2>center(2));
soma2a = mean([T3.soma_center_1(index2a),T3.soma_center_2(index2a)],1);
index2b = (T3.soma_center_1>center(1));
soma2b = mean([T3.soma_center_1(index2b),T3.soma_center_2(index2b)],1);
index2c = (T3.soma_center_1<center(1) & T3.soma_center_2<center(2));
soma2c = mean([T3.soma_center_1(index2c),T3.soma_center_2(index2c)],1);
%%
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
figure;
ax1 = subplot(1,1,1);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [soma2a;soma2b;soma2c];
axon_bias = [real(z),imag(z)];
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
scale1 = 50;
for i = 1:3
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1,...
        soma_center1(i,1)+axon_bias(i,1)*scale1],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1,...
        soma_center1(i,2)+axon_bias(i,2)*scale1],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1),center(2),36,'*','k');
xlim([0,600]);
ylim([0,1200]);
positions = roi.Position;
positions(end+1,:) = positions(1,:);
plot(positions(:,1),positions(:,2),'k','lineWidth',2);
hold on;
plot([center(1),center(1)],[center(2)+200,center(2)-200],'k--','lineWidth',2);
hold on;
plot([center(1),center(1)-200],[center(2),center(2)],'k--','lineWidth',2);