folder = 'E:\matt_simulation';
filename = {'Mirror8_replicated_normed_onequad_seed','NoMirror8_replicated_normed_onequad_seed'};
%%
pwAll_mirror = [];
for i = 1:450
    fname = [filename{1} num2str(i-1,'%03.f')];
    fprintf([fname '\n']);
    load(fullfile(folder,'July20_mirror_detection',fname));
    pwAll(:,5) = pwAll(:,5)+i*80;
    pwAll_mirror = [pwAll_mirror;pwAll];
end
%%
pwAll_nomirror = [];
for i = 1:450
    fname = [filename{2} num2str(i-1,'%03.f')];
    fprintf([fname '\n']);
    load(fullfile(folder,'July20_mirror_detection',fname));
    pwAll(:,5) = pwAll(:,5)+i*80;
    pwAll_nomirror = [pwAll_nomirror;pwAll];
end
%%
% warpingY = 0.75; % affine trasnform scale factor to match illustrtaion
warpingY = 1;
hist_bin = 20;
[unique_spirals_mirror,scolor,low_color_bound,high_color_bound] = density_color_plot(pwAll_mirror,hist_bin);
unique_spirals_mirror_unit = unique_spirals_mirror(:,3)/(20*20*81*150)*35;
[unique_spirals_nomirror,scolor,low_color_bound,high_color_bound] = density_color_plot(pwAll_nomirror,hist_bin);
unique_spirals_nomirror_unit = unique_spirals_nomirror(:,3)/(20*20*81*150)*35;
max_color = max([unique_spirals_mirror_unit; unique_spirals_nomirror_unit]);

h2 = figure;
subplot(1,2,1);
scatter(unique_spirals_nomirror(:,1),unique_spirals_nomirror(:,2)*warpingY,3,unique_spirals_nomirror_unit,'filled');
axis image; 
% axis off;
xlim([0,200]); ylim([0,200]);
title('NoMirror');
caxis([0,0.015]);
colorbar;
colormap(hot);
hold on;
plot([100,100],[0,200],'--w');
hold on;
% plot([0,200],[100,100],'--k');
plot([0,200],[67*2,67*2],'--w');

subplot(1,2,2);
scatter(unique_spirals_mirror(:,1),unique_spirals_mirror(:,2)*warpingY,3,unique_spirals_mirror_unit,'filled');
axis image; 
% axis off;
xlim([0,200]); ylim([0,200]);
title('Mirror');
caxis([0,0.015]);
colorbar;
colormap(hot);
hold on;
plot([100,100],[0,200],'--w');
hold on;
% plot([0,200],[100,100],'--k');
plot([0,200],[67*2,67*2],'--w');
%%
print(h2, 'spiral_detection_mirror_nomirror_simulation-3', '-dpdf', '-bestfit', '-painters');
%%
trialSeq = randperm(450);
trialSeq2 = reshape(trialSeq,[5,90]);
%%
lim1 = 0:10:70;
lim2 = 10:10:80;
spiralN_mirror = zeros(5,8);
ratio_mirror = zeros(5,8);
spiralN_nomirror = zeros(5,8);
ratio_nomirror = zeros(5,8);
for block = 1:5
    trials = trialSeq2(block,:);
    for kk = 1:8
        pwAll_mirror = [];
        for i = 1:90
            itrial = trials(i);
            fname = [filename{1} num2str(itrial-1,'%03.f')];
            fprintf([fname '\n']);
            load(fullfile(folder,'July20_mirror_detection',fname));
            pwAll = pwAll(pwAll(:,5)>lim1(kk) & pwAll(:,5)<lim2(kk),:);
            pwAll(:,5) = pwAll(:,5)+i*80;
            pwAll_mirror = [pwAll_mirror;pwAll];
        end
        pwAll_nomirror = [];
        for i = 1:90
            itrial = trials(i);
            fname = [filename{2} num2str(itrial-1,'%03.f')];
            fprintf([fname '\n']);
            load(fullfile(folder,'July20_mirror_detection',fname));
            pwAll = pwAll(pwAll(:,5)>lim1(kk) & pwAll(:,5)<lim2(kk),:);
            pwAll(:,5) = pwAll(:,5)+i*80;
            pwAll_nomirror = [pwAll_nomirror;pwAll];
        end
        %%
        pwAll_mirror_LD = pwAll_mirror(pwAll_mirror(:,1)<100 & pwAll_mirror(:,2)<132,:);
        pwAll_mirror_RD = pwAll_mirror(pwAll_mirror(:,1)>100 & pwAll_mirror(:,2)<132,:);
        [MatchPoint1,MatchPoint2,DirDiff] = UniqueMatchPoints2(pwAll_mirror_LD,pwAll_mirror_RD);
        sprialAll1_mirror(block,kk) = size(pwAll_mirror_LD,1);
        sprialAll2_mirror(block,kk) = size(pwAll_mirror_RD,1);
        spiralN_mirror(block,kk) = size(MatchPoint1,1);
        ratio_mirror(block,kk) = sum(DirDiff==0)./numel(DirDiff); %match
        %%
        pwAll_nomirror_LD = pwAll_nomirror(pwAll_nomirror(:,1)<100 & pwAll_nomirror(:,2)<132,:);
        pwAll_nomirror_RD = pwAll_nomirror(pwAll_nomirror(:,1)>100 & pwAll_nomirror(:,2)<132,:);
        [MatchPoint1_nomirror,MatchPoint2_nomirror,DirDiff_nomirror] = UniqueMatchPoints2(pwAll_nomirror_LD,pwAll_nomirror_RD);
        sprialAll1_nomirror(block,kk) = size(pwAll_nomirror_LD,1);
        sprialAll2_nomirror(block,kk) = size(pwAll_nomirror_RD,1);
        spiralN_nomirror(block,kk) = size(MatchPoint1_nomirror,1);
        ratio_nomirror(block,kk) = sum(DirDiff_nomirror==0)./numel(DirDiff_nomirror); %match
    end
end
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
color1 = cbrewer2('seq','greys',8);
spiral_rate = spiralN_nomirror./(10*90);
spiral_rate_mean = mean(spiral_rate,1);
spiral_rate_sem = std(spiral_rate,1)./sqrt(5);
ratio_mean = mean(ratio_nomirror,1);
ratio_sem = std(ratio_nomirror,1)./sqrt(5);
subplot(2,2,1);
for block = 1:5
    plot(1:8,spiral_rate,'k');
    hold on;
end
errorbar(1:8,spiral_rate_mean,spiral_rate_sem,'b','lineWidth',2);
ylim([0,0.4]);
xlabel('Iterative time steps');
ylabel('Rotating wave rates');
subplot(2,2,2);
for block = 1:5
    plot(1:8,ratio_nomirror,'k');
    hold on;
end
errorbar(1:8,ratio_mean,ratio_sem,'b','lineWidth',2);
ylim([0.4,1.1]);
xlabel('Iterative time steps');
ylabel('Matching ratio');

spiral_rate = spiralN_mirror./(10*90);
spiral_rate_mean = mean(spiral_rate,1);
ratio_mean = mean(ratio_mirror,1);
spiral_rate_sem = std(spiral_rate,1)./sqrt(5);
ratio_sem = std(ratio_mirror,1)./sqrt(5);
subplot(2,2,3);
for block = 1:5
    plot(1:8,spiral_rate,'k');
    hold on;
end
errorbar(1:8,spiral_rate_mean,spiral_rate_sem,'r','lineWidth',2);
ylim([0,0.4]);
xlabel('Iterative time steps');
ylabel('Rotating wave rates');

subplot(2,2,4);
for block = 1:5
    plot(1:8,ratio_mirror,'k');
    hold on;
end
errorbar(1:8,ratio_mean,ratio_sem,'r','lineWidth',2);
ylim([0.4,1.1]);
xlabel('Iterative time steps');
ylabel('Matching ratio');
print(h1,'matching_rotating_waves_model_single.pdf','-dpdf', '-bestfit', '-painters');
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 200]);
color1 = cbrewer2('seq','greys',8);
spiralAll1a = sprialAll1_mirror./(10*90);
spiralAll_mean = mean(spiralAll1a,1);
spiralAll_sem = std(spiralAll1a,1)./sqrt(5);

spiral_rate = spiralN_nomirror./(10*90);
spiral_rate_mean = mean(spiral_rate,1);
spiral_rate_sem = std(spiral_rate,1)./sqrt(5);

ratio_mean = mean(ratio_nomirror,1);
ratio_sem = std(ratio_nomirror,1)./sqrt(5);

subplot(1,3,1);
errorbar(1:8,spiralAll_mean,spiralAll_sem,'k','lineWidth',1);
subplot(1,3,2);
errorbar(1:8,spiral_rate_mean,spiral_rate_sem,'k','lineWidth',1);
subplot(1,3,3);
errorbar(1:8,ratio_mean,ratio_sem,'k','lineWidth',1);

spiralAll1a = sprialAll1_nomirror./(10*90);
spiralAll_mean = mean(spiralAll1a,1);
spiralAll_sem = std(spiralAll1a,1)./sqrt(5);

spiral_rate = spiralN_mirror./(10*90);
spiral_rate_mean = mean(spiral_rate,1);
ratio_mean = mean(ratio_mirror,1);
spiral_rate_sem = std(spiral_rate,1)./sqrt(5);
ratio_sem = std(ratio_mirror,1)./sqrt(5);
subplot(1,3,1);
hold on;
errorbar(1:8,spiralAll_mean,spiralAll_sem,'r','lineWidth',1);
ylim([0,1.5]);
yticks([0:0.5:1.5]);
yticklabels([0:0.5:1.5]);
xlabel('Iterative time steps');
ylabel('Rotating wave rates');
box off;

subplot(1,3,2);
hold on;
errorbar(1:8,spiral_rate_mean,spiral_rate_sem,'r','lineWidth',1);
% ylim([0,0.4]);
% yticks([0:0.1:0.4]);
% yticklabels([0:0.1:0.4]);
ylim([0,1.5]);
yticks([0:0.5:1.5]);
yticklabels([0:0.5:1.5]);
xlabel('Iterative time steps');
ylabel('Rotating wave rates');
legend({'NoMirror','Mirror'});
box off;

subplot(1,3,3);
hold on;
errorbar(1:8,ratio_mean,ratio_sem,'r','lineWidth',1);
ylim([0.4,1.1]);
yticks([0.4:0.2:1.0]);
yticklabels([0.4:0.2:1.0]);
xlabel('Iterative time steps');
ylabel('Matching ratio');
box off;
%%
print(h1,'matching_rotating_waves_model_group.pdf','-dpdf', '-bestfit', '-painters');