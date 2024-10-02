load(fullfile(data_folder,'spirals\spirals_density',...
    'spiralDensityLinePerSession.mat'));                                   % desnity line across session
%%
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_fft_bump');
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(td,'yyyymmdd');
    subfolder = [mn '_' tdb '_' num2str(en)];
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,[fname '_histogram.mat']));
    count_sample_alpha_all(:,kk) = count_sample_alpha;
    count_sample_nonalpha_all(:,kk) = count_sample_nonalpha;
end
%%
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
%%
figure;
for kk = 1:15
    subplot(3,5,kk);
    plot(points_line,count_sample_alpha_all(:,kk),'r');
    hold on;
    plot(points_line,count_sample(:,kk),'k');
end
