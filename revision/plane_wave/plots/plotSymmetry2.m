function hr1b = plotSymmetry2(data_folder,save_folder)
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
ax3 = subplot(2,4,1); 
bar(1:3,[ratio_mean(1,3:4),ratioM_mean(1,2)],1,'faceColor',[0.5,0.5,0.5]);
hold on;
er = errorbar(1:3,[ratio_mean(1,3:4),ratioM_mean(1,2)],[ratio_sem(1,3:4),ratioM_sem(1,2)]);
er.Color = 'k';
er.LineStyle = 'None';
ylim([0,0.4]);
xticks(1:3);
xticklabels({'SSp-L','SSp-R','SSp-LR'});

ax3 = subplot(2,4,5); 
bar(1:3,[ratio_mean(1,[2,4]),ratioM_mean(1,1)],1,'faceColor',[0.5,0.5,0.5]);
hold on;
er = errorbar(1:3,[ratio_mean(1,[2,4]),ratioM_mean(1,1)],[ratio_sem(1,[2,4]),ratioM_sem(1,1)]);
er.Color = 'k';
er.LineStyle = 'None';
ylim([0,0.4]);
xticks(1:3);
xticklabels({'MO-R','SSp-R','MO-SSp-R'});
%
plots = [3,4,2,4];
for k = 1:4
    if ismember(k,[3,4])
        k1 = k+3;
    else
        k1 = k+1;
    end
    ax3 = subplot(2,4,k1); 
    amp_ssp_r = squeeze(amp_all(:,:,4));
    angle_temp = squeeze(angle_all(:,:,plots(k)));
    amp_temp = squeeze(amp_all(:,:,plots(k)));
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
    ylim([0,0.4]);
end

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
subplot(2,4,4);
scatter(angle1,angle2,4,'k','filled');
xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylim([-pi,pi]);
yticks([-pi:pi/2:pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    
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
subplot(2,4,8);
scatter(angle1,angle2,4,'k','filled');
xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ylim([-pi,pi]);
yticks([-pi:pi/2:pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
%%
print(hr1b, fullfile(save_folder,'FigR1b_wave_symmetry.pdf'),...
    '-dpdf', '-bestfit', '-painters');