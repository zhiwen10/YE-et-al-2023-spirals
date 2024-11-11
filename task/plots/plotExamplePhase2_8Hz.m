function hs14efg = plotExamplePhase2_8Hz(data_folder,save_folder)
%%
freq_high = [2,8];
freq_low = [0.05,2];
clear wf_all phase_all amp_all contrast_all
mn ='ZYE_0091';
trialWin = [-4,4];
load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
load(fullfile(data_folder,'task','task_outcome',[mn '_task_freq_to' num2str(freq_high(2)) 'Hz']));
amp_high = amp_all;
load(fullfile(data_folder,'task','task_outcome',[mn '_task_freq_to' num2str(freq_low(2)) 'Hz']),'amp_all');
amp_high1 = mean(squeeze(amp_high(123:141,1,:)),1)';
p_high1 = prctile(amp_high1,75);
contrasts = [0.06, 0.125, 0.25, 0.5, 1,0,0,0,0,0;  ...
    0, 0, 0, 0, 0, 0.06,0.125,0.25,0.5,1];
t1 = [trialWin(1):1/35:trialWin(2)]';
stimOn = (numel(t1)-1)/2+1;
wf_all = wf_all*100;
%%
xlim1 = 2;
hs14efg = figure('Renderer', 'painters', 'Position', [50  50 850 550]);
for i = 1:5
    indx = (T_all.left_contrast == contrasts(1,i) & ...
        T_all.right_contrast == contrasts(2,i) & T_all.label == "correct");
    wf_temp = wf_all(:,1:2,indx);
    wf_mean = squeeze(mean(wf_temp,3));
    %%
    indx1 = (T_all.left_contrast == contrasts(1,i) & ...
        T_all.right_contrast == contrasts(2,i) & T_all.label == "correct" & amp_high1>p_high1);
    wf_temp1 = wf_all(:,1:2,indx1);
    %%
    subplot(5,10,(i-1)*2+1);
    plot(t1,squeeze(wf_mean(:,1)),'k'); % left hemisphere
    xline(t1(stimOn),'k--');
    xlim([-xlim1,xlim1]);
    ylim([-3,3]);
    
    subplot(5,10,(i-1)*2+2);
    plot(t1,squeeze(wf_mean(:,2)),'r'); % right hemisphere
    xline(t1(stimOn),'k--');
    xlim([-xlim1,xlim1]);
    ylim([-3,3]);
    %%       
    subplot(5,10,[(i-1)*2+11,(i-1)*2+21]);
    for k = 1:20
        plot(t1,squeeze(wf_temp1(:,1,k))+k*6,'k'); % left hemisphere
        hold on;
    end
    xline(t1(stimOn),'k--');
    xlim([-xlim1,xlim1]);
    ylim([-20,140]);
    set(gca,'YTickLabel',[]);
    if i ==1
        hold on;
        plot([0, 0],[0,8],'r');
    end
    
    subplot(5,10,[(i-1)*2+12,(i-1)*2+22]);
    for k = 1:20
        plot(t1,squeeze(wf_temp1(:,2,k))+k*6,'r'); % right hemisphere
        hold on;
    end
    xline(t1(stimOn),'k--');
    xlim([-xlim1,xlim1]);
    ylim([-20,140]);
    set(gca,'YTickLabel',[]);
    %%
    phase_temp = squeeze(phase_all(141,1:2,indx));
    subplot(5,10,[(i-1)*2+31:(i-1)*2+32]);
    scatter(phase_temp(1,:),phase_temp(2,:),6,'k','filled'); % left hemisphere
    axis equal;
    xlim([-pi,pi]);
    ylim([-pi,pi]);
    xticks([-pi,-pi/2,0,pi/2,pi]);
    xticklabels({'-pi','-pi/2','0','pi/2','pi'});
    yticks([-pi,-pi/2,0,pi/2,pi]);
    yticklabels({'-pi','-pi/2','0','pi/2','pi'});
    %%
    subplot(5,10,[(i-1)*2+41:(i-1)*2+42]);
    edges = [-pi:pi/4:pi];
    histogram(phase_temp(1,:),edges,'EdgeColor','k','FaceColor',[0.5,0.5,0.5]);
    xlim([-pi,pi]);
    ylim([0,150]);
    xticks([-pi,-pi/2,0,pi/2,pi]);
    xticklabels({'-pi','-pi/2','0','pi/2','pi'});
end
%%
print(hs14efg,fullfile(save_folder,'FigS14efg_phase_histogram_2_8Hz'),'-dpdf', '-bestfit', '-painters');
