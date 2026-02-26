h1 = figure('Renderer', 'painters', 'Position', [100 100 700 400]);
subplot(1,2,1);
load('cutting_counts2.mat');
spiral_density_control2a = spiral_density_control2;
for i = 1:3
    scatter(ones(size(spiral_density_control2a))*(2*i-1),spiral_density_control2a(:,i),8,'k');
    hold on;
end
for i = 1:3
    scatter(ones(size(spiral_density_control2a))*2*i,spiral_density_cutting2(:,i),8,'r');
    data_temp1 = spiral_density_control2a(:,i);
    data_temp2 = spiral_density_cutting2(:,i);
    plot([ones(size(data_temp1))*(2*i-1),ones(size(data_temp2))*(2*i)]',...
        [data_temp1,data_temp2]','color','k');
end
xticks([1.5:2:7.5]);
xticklabels({'MOs','MOp','SSp'});
yticks([0:5:15]);
ylim([0,15]);
ylabel('Total spirals/s');
[hh1,pp1] = ttest(spiral_density_control2a(:,1),spiral_density_cutting2(:,1));
[hh2,pp2] = ttest(spiral_density_control2a(:,2),spiral_density_cutting2(:,2));
[hh3,pp3] = ttest(spiral_density_control2a(:,3),spiral_density_cutting2(:,3));

control2a_mean = mean(spiral_density_control2a,1);
control2a_sem = std(spiral_density_control2a,1)./sqrt(5);
cutting2a_mean = mean(spiral_density_cutting2,1);
cutting2a_sem = std(spiral_density_cutting2,1)./sqrt(5);

subplot(1,2,2);
load('control_counts2.mat');
% spiral_density_control2(4,:) = [];
% spiral_density_craniotomy2(4,:) = [];

for i = 1:3
    scatter(ones(size(spiral_density_control2))*(2*i-1),spiral_density_control2(:,i),8,'k');
    hold on;
end
for i = 1:3
    scatter(ones(size(spiral_density_control2))*2*i,spiral_density_craniotomy2(:,i),8,'r');
    data_temp1 = spiral_density_control2(:,i);
    data_temp2 = spiral_density_craniotomy2(:,i);
    plot([ones(size(data_temp1))*(2*i-1),ones(size(data_temp2))*(2*i)]',...
        [data_temp1,data_temp2]','color','k');
end
xticks([1.5:2:7.5]);
xticklabels({'MOs','MOp','SSp'});
yticks([0:5:15]);
ylim([0,15]);
ylabel('Total spirals/s');

[hh1b,pp1b] = ttest(spiral_density_control2(:,1),spiral_density_craniotomy2(:,1));
[hh2b,pp2b] = ttest(spiral_density_control2(:,2),spiral_density_craniotomy2(:,2));
[hh3b,pp3b] = ttest(spiral_density_control2(:,3),spiral_density_craniotomy2(:,3));

control2_mean = mean(spiral_density_control2,1);
control2_sem = std(spiral_density_control2,1)./sqrt(10);
craniotomy2_mean = mean(spiral_density_craniotomy2,1);
craniotomy2_sem = std(spiral_density_craniotomy2,1)./sqrt(10);
%%
print(h1, 'Cutting_control_summary2', '-dpdf', '-bestfit', '-painters');  
%%
ratio_cutting =  (spiral_density_control2a-spiral_density_cutting2)./spiral_density_control2a;
ratio_cutting_median = median(ratio_cutting,1);
ratio_cutting_sem = std(ratio_cutting,[],1)./sqrt(size(ratio_cutting,1));
ratio_control =  (spiral_density_control2-spiral_density_craniotomy2)./spiral_density_control2;
ratio_control_median = median(ratio_control,1);
ratio_control_sem = std(ratio_control,[],1)./sqrt(size(ratio_control,1));
%%
ratio_cutting_mean2 = (mean(spiral_density_control2a,1)-mean(spiral_density_cutting2,1))./mean(spiral_density_control2a,1);
ratio_control_mean2 =  (mean(spiral_density_control2,1)-mean(spiral_density_craniotomy2,1))./mean(spiral_density_control2,1);
