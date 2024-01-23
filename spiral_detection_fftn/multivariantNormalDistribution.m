mu = [3 3];
% Sigma = [1 1.5; 1.5 3];
Sigma = eye(2);
rng('default')  % For reproducibility
R = mvnrnd(mu,Sigma,2000);
figure;
plot(R(:,1),R(:,2),'+');
%%
[f,xi] = ksdensity([R(:,1),R(:,2)],[R(:,1),R(:,2)]); 
figure
scatter(R(:,1), R(:,2), [], f);
%%
F = scatteredInterpolant(R(:,1), R(:,2),f);
a = -1:0.2:7;
[xq,yq] = meshgrid(a);
vq1 = F(xq,yq)*2000;
figure;surf(xq,yq,vq1)
%% integrate to get area under the curve, result is close to 1
I1 = trapz(a ,trapz(a ,vq1,2));
%%
figure;
scatter_kde(R(:,1),R(:,2))
colorbar
%%
% Generate data
x = normrnd(10,1,1000,1);
y = normrnd(10,1,1000,1);
% Plot data using probability density estimate as function
figure(1); 
h1 = scatter_kde(x, y, 'filled', 'MarkerSize', 100);
% Add Color bar
cb = colorbar();
cb.Label.String = 'Probability density estimate';
%%
rng('default')  % For reproducibility
xx = [0+.5*rand(20,1) 5+2.5*rand(20,1);
            .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
%%
gridx1 = -0.25:.05:1.25;
gridx2 = 0:.1:15;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
X = [0+.5*rand(20,1) 5+2.5*rand(20,1);
  .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
ksdensity(X,xi);