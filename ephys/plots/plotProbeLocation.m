function h4f = plotProbeLocation(T,data_folder,save_folder)
mainfolder  = fullfile(data_folder,'ephys','probe_location');
[atlas, metaAVGT] = nrrdread(fullfile(data_folder,'tables',...
    'annotation_50.nrrd'));
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
% area = {'CORTEX','THAL','STR','MB','PONS'};
% color_id = [2,3,1,4,9]; % blue, green, red, purple,grey in sequence of color1
area = {'CORTEX','THAL','STR','MB'};
color_id = [2,3,1,4]; % blue, green, red, purple,grey in sequence of color1
color1 = cbrewer2('qual','Set1',9);
%%
% 3 subareas
root1 = {'/997/'};
maskPath{1}  = '/997/8/343/1129/549/'; % Th
maskPath{2}  = '/997/8/567/623/477/'; % STR
maskPath{3}  = '/997/8/343/313/'; % midbrain
maskPath{4}  = "/997/8/567/688/695/315/"; % isocortex
maskPath{5}  = "/997/8/567/688/695/698/"; % olfactory

atlas1{1} = double(squeeze(atlas(80,:,:)))';
atlas1{2} = double(squeeze(atlas(:,:,80)));
atlas1{3} = double(squeeze(atlas(:,120,:)));

order_2d(:,1) = [3,1];
order_2d(:,2) = [3,2];
order_2d(:,3) = [1,2];
axis_lim(:,1) = [264,264];
axis_lim(:,2) = [264,264];
axis_lim(:,3) = [264,264];
%% plot
h4f = figure('Renderer', 'painters', 'Position', [100 100 300 600]);
for iplot = 1:3
    subplot(3,1,iplot);
    scale3 = 1;
    plotOutline_fill(root1,st,atlas1{iplot},[],scale3,'w');
    plotOutline_fill(maskPath(1),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
    plotOutline_fill(maskPath(2),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
    plotOutline_fill(maskPath(3),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
    plotOutline_fill(maskPath(4:5),st,atlas1{iplot},[],scale3,'w');
    axis image; axis off;
    set(gca, 'YDir','reverse');
    xlim([-20,axis_lim(1,iplot)]+20);
    ylim([-20,axis_lim(2,iplot)]+20);
    for current_area = [2,3,4]
        %%
        indx = find(contains(T.Area,area{current_area}));     
        current_T = T(indx,:);
        count = 1;
        atlas_convert = permute(atlas,[2,1,3]);
        clear points pointColor lineColor
        for kk = 1:size(current_T,1)
            ops = get_session_info2(current_T,kk,data_folder);
            fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
            probefolder_full = fullfile(mainfolder,fname);
            pointDir = dir(fullfile(probefolder_full, '*.csv'));
            for j = 1:numel(pointDir)
                pointcsv = fullfile(probefolder_full,pointDir(j).name);
                pointRaw = csvread(pointcsv);
                pointRaw = round(pointRaw/2);
                pointRaw = pointRaw(:,[2,1,3]);
                if current_area==4 & kk>=5
                    pointRaw(:,1) = 228-pointRaw(:,1);
                end
                [m,p,s] = best_fit_line(pointRaw(:,1), pointRaw(:,2), pointRaw(:,3));
                % ensure proper orientation: want 0 at the top of the brain
                % and positive distance goes down into the brain
                if p(2)<0
                    p = -p;
                end
                %%
                % determine "origin" at top of brain
                % step upwards along tract direction until tip of brain / past cortex
                ann = 50;
                in_brain = true;
                vector = m-pointRaw(1,:);
                vector_length_old = vecnorm(vector,2,2);           
                while in_brain
                    m = m-p/5;
                    vector1 = m-pointRaw(1,:);
                    vector_length_new = vecnorm(vector1,2,2);
                    if vector_length_new>vector_length_old
                        in_brain = false;
                    end
                    vector_length_old = vector_length_new;    
                end
                [depth_max, tip_index_max] = max(pointRaw(:,2));
                [depth_min, tip_index_min] = min(pointRaw(:,2));
                reference_probe_length_tip = sqrt(sum(...
                    (pointRaw(tip_index_max,:) - pointRaw(tip_index_min,:)).^2));

                % display the scaling
                % plot line the length of the probe in reference space
                probe_length_histo = round(reference_probe_length_tip);
                hold on;
                % plot line the length of the entire probe in reference space
                plot(m(order_2d(1,iplot))+p(order_2d(1,iplot))*...
                    [1 probe_length_histo],...
                    m(order_2d(2,iplot))+p(order_2d(2,iplot))*...
                    [1 probe_length_histo], ...
                    'Color', color1(color_id(current_area),:), ...
                    'LineWidth', 1);
                hold on;
            end
        end
    end
end
%%
print(h4f, fullfile(save_folder,'Fig4f_probe_projection_2d'),...
    '-dpdf', '-bestfit', '-painters');
