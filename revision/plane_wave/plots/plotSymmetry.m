load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all.mat'));
%%
vxy_all = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
angle_all = angle(vxy_all);
amp_all = abs(vxy_all);
%%
figure;
for k = 1:4
    ax3 = subplot(1,4,k); 
    amp_ssp_r = squeeze(amp_all(:,:,4));
    angle_temp = squeeze(angle_all(:,:,k));
    amp_temp = squeeze(amp_all(:,:,k));
    edges = [-pi:pi/6:pi];
    for i = 1:15
        angle_temp1 = angle_temp(i,:);
        amp_temp1 = amp_temp(i,:);
        amp_ssp_r1 = amp_ssp_r(i,:);
        angle_temp2 = angle_temp1(amp_temp1>=0.6 & amp_ssp_r1>=0.6);
        [N(i,:),edges] =histcounts(angle_temp2,edges);
        N1(i,:) = N(i,:)/size(angle_temp2,2);
    end
    N_mean = mean(N1,1);
    N_sem = std(N1,1)./sqrt(15);

    N_mean(end+1) = N_mean(1);
    N_sem(end+1) = N_sem(1);
    bar(edges,N_mean,1,'faceColor',[0.5,0.5,0.5]);
    hold on;
    er = errorbar(edges,N_mean,N_sem);
    % er = errorbar(edges+pi/12,N_mean,N_sem);
    er.Color = 'k';
    er.LineStyle = 'None';
    % xlim([-pi,pi]);
    xticks([-pi:pi/2:pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    ylim([0,0.35]);
end
%%
for i = 1:15 
    clear angle_temp1 index angle_temp2 angle2
    angle_temp1 = squeeze(angle_all(i,:,4));
    amp_temp1 = squeeze(amp_all(i,:,4));
    angle_temp2 = squeeze(angle_all(i,:,2));
    amp_temp2 = squeeze(amp_all(i,:,2));
    % index = (amp_temp1>=0.6 & amp_temp2>=0.6 & angle_temp1>=0 & angle_temp1<=pi);    
    index = (amp_temp1>=0.6 & amp_temp2>=0.6);  
    angle1 = angle_temp1(index);
    angle2 = angle_temp2(index);
    [N2(i,:),edges] =histcounts(angle2,edges);
    N2a(i,:) = N2(i,:)/size(angle2,2);
end
figure;
N_mean = mean(N2a,1);
N_sem = std(N2a,1)./sqrt(15);

N_mean(end+1) = N_mean(1);
N_sem(end+1) = N_sem(1);
bar(edges,N_mean,1,'faceColor',[0.5,0.5,0.5]);
hold on;
er = errorbar(edges,N_mean,N_sem);
% er = errorbar(edges+pi/12,N_mean,N_sem);
er.Color = 'k';
er.LineStyle = 'None';
% xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylim([0,0.5]);
%%
angle_all1 = squeeze(angle_all(:,:,4));
angle_all1 = angle_all1(:);
angle_all2 = squeeze(angle_all(:,:,2));
angle_all2 = angle_all2(:);
amp_all1 = squeeze(amp_all(:,:,4));
amp_all1 = amp_all1(:);
amp_all2 = squeeze(amp_all(:,:,2));
amp_all2 = amp_all2(:);
index = (amp_all1>=0.6 & amp_all2>=0.6);  
angle1 = angle_all1(index);
angle2 = angle_all2(index);
figure;
scatter(angle1,angle2);
%%
angle_all1 = squeeze(angle_all(:,:,4));
angle_all1 = angle_all1(:);
angle_all2 = squeeze(angle_all(:,:,3));
angle_all2 = angle_all2(:);
amp_all1 = squeeze(amp_all(:,:,4));
amp_all1 = amp_all1(:);
amp_all2 = squeeze(amp_all(:,:,3));
amp_all2 = amp_all2(:);
index = (amp_all1>=0.6 & amp_all2>=0.6);  
angle1 = angle_all1(index);
angle2 = angle_all2(index);
figure;
scatter(angle1,angle2);


