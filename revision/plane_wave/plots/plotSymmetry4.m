function hs8o = plotSymmetry4(data_folder,save_folder)
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all.mat'));
vxy_all = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
angle_all = angle(vxy_all);
amp_all = abs(vxy_all);
%%
threthold = 0.6;
for k = 1:4
    angle_all1 = squeeze(angle_all(:,:,k));
    amp_all1 = squeeze(amp_all(:,:,k));
    for i = 1:15
        amp_temp = amp_all1(i,:);
        ratio(i,k) = sum(amp_temp>= threthold)./numel(amp_temp);
    end
end
ratio_mean = mean(ratio,1);
ratio_sem = std(ratio,1)./sqrt(15);
%%
for k = 1:2
    angle_all1 = squeeze(angle_all(:,:,k+1));
    amp_all1 = squeeze(amp_all(:,:,k+1));
    amp_all2 = squeeze(amp_all(:,:,4));
    for i = 1:15
        amp_temp = amp_all1(i,:);
        amp_temp2 = amp_all2(i,:);
        ratioM(i,k) = sum(amp_temp>threthold & amp_temp2 >threthold)./numel(amp_temp);
    end
end
ratioM_mean = mean(ratioM,1);
ratioM_sem = std(ratioM,1)./sqrt(15);
%%
hs8o = figure('Renderer', 'painters', 'Position', [50 50 700 300]);
subplot(1,2,1);
angle_all1 = squeeze(angle_all(:,:,4));
angle_all1 = angle_all1(:);
angle_all2 = squeeze(angle_all(:,:,3));
angle_all2 = angle_all2(:);
amp_all1 = squeeze(amp_all(:,:,4));
amp_all1 = amp_all1(:);
amp_all2 = squeeze(amp_all(:,:,3));
amp_all2 = amp_all2(:);
index = (amp_all1>=threthold & amp_all2>=threthold);  
angle1 = angle_all1(index);
angle2 = angle_all2(index);
scatter(angle1,angle2,4,'k','filled');
xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylim([-pi,pi]);
yticks([-pi:pi/2:pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    
subplot(1,2,2);
angle_all1 = squeeze(angle_all(:,:,4));
angle_all1 = angle_all1(:);
angle_all2 = squeeze(angle_all(:,:,2));
angle_all2 = angle_all2(:);
amp_all1 = squeeze(amp_all(:,:,4));
amp_all1 = amp_all1(:);
amp_all2 = squeeze(amp_all(:,:,2));
amp_all2 = amp_all2(:);
index = (amp_all1>=threthold & amp_all2>=threthold);  
angle1 = angle_all1(index);
angle2 = angle_all2(index);
scatter(angle1,angle2,4,'k','filled');
xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylim([-pi,pi]);
yticks([-pi:pi/2:pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%%
print(hs8o, fullfile(save_folder,'FigS8o_wave_symmetry.pdf'),...
    '-dpdf', '-bestfit', '-painters');