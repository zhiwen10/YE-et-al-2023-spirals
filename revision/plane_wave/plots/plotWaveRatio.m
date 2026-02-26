function hr1b = plotWaveRatio(data_folder,save_folder)
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all.mat'));
vxy_all = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
angle_all = angle(vxy_all);
amp_all = abs(vxy_all);
%%
for k = 1:4
    angle_all1 = squeeze(angle_all(:,:,k));
    amp_all1 = squeeze(amp_all(:,:,k));
    for i = 1:15
        amp_temp = amp_all1(i,:);
        ratio(i,k) = sum(amp_temp>0.6)./numel(amp_temp);
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
        ratioM(i,k) = sum(amp_temp>0.6 & amp_temp2 >0.6)./numel(amp_temp);
    end
end
ratioM_mean = mean(ratioM,1);
ratioM_sem = std(ratioM,1)./sqrt(15);
%%
hr1b = figure('Renderer', 'painters', 'Position', [50 50 700 350]);
ax3 = subplot(2,2,1); 
bar(1:4,ratio_mean,1,'faceColor',[0.5,0.5,0.5]);
hold on;
er = errorbar(1:4,ratio_mean,ratio_sem);
er.Color = 'k';
er.LineStyle = 'None';
ylim([0,0.4]);
xticks(1:4);
xticklabels({'MO-L','MO-R','SSp-L','SSp-R'});
%%
ax3 = subplot(2,2,2); 
bar(1:3,[ratio_mean(1,[2,4]),ratioM_mean(1,1)],1,'faceColor',[0.5,0.5,0.5]);
hold on;
er = errorbar(1:3,[ratio_mean(1,[2,4]),ratioM_mean(1,1)],[ratio_sem(1,[2,4]),ratioM_sem(1,1)]);
er.Color = 'k';
er.LineStyle = 'None';
ylim([0,0.4]);
xticks(1:3);
xticklabels({'MO-L','MO-R','SSp-L','SSp-R'});