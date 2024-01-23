h = figure('Renderer', 'painters', 'Position', [100 100 600 700]);
subplot(1,1,1);
load('example_12_cells.mat');
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for icell = 1
    axon_current_all = axon_current_all1{icell};
    axon_current = axon_current1{icell};
    soma_current_2d = soma_current_2d1{icell};
    
    axon_current_2d = double(axon_current(:,[3,1]));
    
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);    
    soma_center1 = double(soma_current_2d(:,[3,1]));
    sample_n = size(axon_current_2d,1);
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color2(angle1,:),'MarkerSize',1); 
    hold on;
    scatter(axon_current(:,3),axon_current(:,1),4,'k','filled');
    hold on;
%     plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
%         [soma_center1(2)-U(2,1)*scale1*polarity3(icell),soma_center1(2)+U(2,1)*scale1*polarity3(icell)],'color','k','LineWidth',1);
    plot([soma_center1(1)-U(1,1)*S_diag3(1)/sqrt(sample_n),soma_center1(1)+U(1,1)*S_diag3(1)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,1)*S_diag3(1)/sqrt(sample_n),soma_center1(2)+U(2,1)*S_diag3(1)/sqrt(sample_n)],'color','r','LineWidth',1);
    hold on;
    plot([soma_center1(1)-U(1,2)*S_diag3(2)/sqrt(sample_n),soma_center1(1)+U(1,2)*S_diag3(2)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,2)*S_diag3(2)/sqrt(sample_n),soma_center1(2)+U(2,2)*S_diag3(2)/sqrt(sample_n)],'color','r','LineWidth',1);
    % axis off; 
    axis image;
    xlim([0,600]);
    ylim([0,1200]);
end
hold on;
center = [244,542]; % SSp-un
scatter(244,542,36,'x','r');
hold on;
plot([soma_center1(1),244],[soma_center1(2),542],'--','color','k');
%%
% mu = [2 3];
% Sigma = [1 1.5; 1.5 3];
% rng('default')  % For reproducibility
% R = mvnrnd(mu,Sigma,100);
%%
% r1 = normrnd(0,1,100,1);
% r2 = normrnd(0,5,100,1);
% R = [r1 r2];
% mean_R = mean(R,1);
% R1 = R-mean_R;
% [U1,S1,V1] = svd(R1,"econ");
% S_diag1 = diag(S1);
% %
% % std = singular_value/sqrt(sample_size)
% figure
% plot(R(:,1),R(:,2),'+');
% hold on;
% plot([mean_R(1)-V1(1,1)*S_diag1(1),mean_R(1)],[mean_R(2)-V1(2,1)*S_diag1(1),mean_R(2)],'k');
% hold on;
% plot([mean_R(1)-V1(1,2)*S_diag1(2),mean_R(1)],[mean_R(2)-V1(2,2)*S_diag1(2),mean_R(2)],'k');
% axis image