function [edges,ratio_all_mean,ratio_all_sem,ratio_all] = ratio_amp(T1,amp_folder,spiral_left_match_all,spiral_right_match_all)
for kk = 1:size(T1,1)
    clear spiral_left_temp spiral_right_temp spiral_all_temp
    ops = get_session_info(T1,kk);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    serverRoot = expPath(ops.mn, ops.td, ops.en);   
    load(fullfile(amp_folder,[fname '_amp.mat']));
    t = readNPY(fullfile(serverRoot, 'blue','svdTemporalComponents.timestamps.npy'));
    spiral_left_temp = spiral_left_match_all{kk};
    spiral_left_temp(:,12) = traceAmp(1,spiral_left_temp(:,5))';
    spiral_right_temp = spiral_right_match_all{kk};
    spiral_right_temp(:,12) = traceAmp(1,spiral_right_temp(:,5))';    
    spiral_all_temp = [spiral_left_temp;spiral_right_temp];
    
    edges = [0:0.0025:0.025];
    % edges = [0:0.00125:0.025];
    for i = 1:numel(edges)-1
        a = find(spiral_all_temp(:,12)>=edges(i) & spiral_all_temp(:,12)<edges(i+1));
        spiral_all_temp1 = spiral_all_temp(a,:);
        if size(spiral_all_temp1,1)>10
            ratio_all(kk,i) = sum(spiral_all_temp1(:,11))/size(spiral_all_temp1,1);
        else
            ratio_all(kk,i) = nan;
        end     
    end
end
ratio_all_mean = mean(ratio_all,1,"omitnan");
ratio_all_std = std(ratio_all,[],1,"omitnan");
count_all = sum(not(isnan(ratio_all)),1);
ratio_all_sem = ratio_all_std./sqrt(count_all);
