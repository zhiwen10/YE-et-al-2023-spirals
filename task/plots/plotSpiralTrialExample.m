function hs15gh = plotSpiralTrialExample(data_folder, save_folder)
%% load atlas brain horizontal projection and outline 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
%%
mn = 'ZYE_0085';
session = 6;
load(fullfile(data_folder,'task','spirals_half','ZYE_0085_half_spirals4.mat'));
load(fullfile(data_folder,'task','spirals','ZYE_0085_spirals_task_sort.mat'),'spiral_all');
load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome.mat']));
load(fullfile(data_folder,'task','task_outcome',[mn '_task_trial_ID.mat']));
%%
index_session = (trial_all.session == session);
spiral_session = spiral_all(index_session,:);
trial_session = trial_all(index_session,:);
%%
T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
T1 = T_session(T_session.label == "task",:);
T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
mn = char(T1.MouseID{session});
tda = T1.date(session);  
en = T1.folder(session); 
block_en = T1.block_folder(session); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
expDir = dir(fullfile(session_root,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
% get photodiode time
win = [0,3000];
allPD2 = getPhotodiodeTime(session_root,win);
%%
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);                    % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1)';
allPD2 = allPD2(1:ntrial);
allPD2 = allPD2(norepeat_indx);
%% load atlas registration
downscale = 16;
load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(dV(1:50,:));
%% get time index for all events
clear indx
for i = 1:numel(allPD2)
    ta  = t-allPD2(i);
    indx2(i,1) = find(ta>0, 1,'first');
    t2(i,1) = t(indx2(i,1));
end
%%
lineColor = 'k';
downscale = 16;
hemi = [];
BW3 = BW(1:downscale:end,1:downscale:end);
scale3 = 5/16;
lineColor = 'k'; lineColor1 = 'w';
%%
indx = [576,623,608,54,65,136]; % 3 spirla example, 3 no spiral example 
cmax = 0.03; cmin = -0.03;
hs15gh =figure('Renderer', 'painters', 'Position', [50 50 950 400]);
for kk = 1:numel(indx)
    current_trial = indx(kk);
    indx_t = indx2(current_trial);
    %%
    Ut1a = double(reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3)));
    trace = Ut1a*dV1(:,indx_t-4:indx_t+32);
    trace2d = reshape(trace,size(Ut1,1),size(Ut1,2),size(trace,2))./mimgt;
    %%
    if kk<=3
        row = kk+1;
    else
        row = kk+2;
    end
    for i = 1:18
        ax4 = subplottight(10,19,(row-1)*19+i);
        im1 = imagesc(squeeze(trace2d(:,:,i)));
        hold on;
        plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
        axis image; axis off;
        set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
        colormap(ax4,'parula');
        caxis([cmin,cmax]);
    end
    
    cb4 = subplottight(10,19,19);
    im1 = imagesc(squeeze(trace2d(:,:,19)),'visible','off');
    colormap(cb4,'parula');
    caxis([cmin,cmax]);
    axis off;
    cb = colorbar;
    a =  cb.Position; %gets the positon and size of the color bar
    set(cb,'Position',[a(1)-0.02 a(2) a(3) a(4)/8])% To change size
    
    if kk == 1
        dim = [.5 .9 .1 .1];
        str = 'Spiral';
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    elseif kk == 4
        dim = [.5 .5 .1 .1];
        str = 'No spiral';
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end
%%
print(hs15gh,fullfile(save_folder,'FigS15gh_example_spiral_nospiral')...
    ,'-dpdf', '-bestfit', '-painters');
