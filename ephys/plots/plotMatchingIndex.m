function [h4h,h4i] = plotMatchingIndex(T,data_folder,save_folder)
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
rfolder = fullfile(data_folder,'ephys','flow_var');
%%
color1 = {'g','r','m'};
sig_ratio = cell(3,1);
flow_var_all = cell(3,1);
count2 = 1;
for current_area = [1,2,4]
    %%
    indx = find(contains(T.Area,area{current_area}));
    current_T = T(indx,:);
    list = 1:size(current_T,1); 
    pa_all = []; ha_all = []; pp_all = []; flow_var_max = [];
    count1 = 1;
    for kk = list
        h = []; p = [];
        ops = get_session_info2(current_T,kk,data_folder);
        fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
        ffname = fullfile(rfolder,fname);
        area_folder = area{current_area};
        load(fullfile(rfolder,[fname '.mat']));
        flow_var1 = 1-flow_var;
        nan_index = isnan(flow_var1(1,:));
        flow_var1(:,nan_index) = [];
        %% unbalanced two-way-anova
        edges1 = edges(1,2:end);
        edges1(nan_index) = [];
        edges_matrix = repmat(edges1,size(flow_var1,1),1);
        permute_id_matrix = ones(size(flow_var1));
        permute_id_matrix(1,:) = 0;        
        flow_var3 = flow_var1(:);
        edges_matrix2 = edges_matrix(:);
        permute_id_matrix2 = permute_id_matrix(:);
        pp = anovan(flow_var3,{edges_matrix2 permute_id_matrix2},...
            'model',2,'varnames',{'amp','permute'},'display','off');
        pp_all = [pp_all;pp'];
        indx1 = not(isnan(flow_var(1,:)));
        for ii = 1:size(flow_var1,2)
            [h(ii),p(ii)] = ttest2(flow_var1(1,ii),flow_var1(2:end,ii),"Tail","right");
        end
        ha = nan(1,20); ha(1,indx1) = h;
        ha_all = [ha_all;ha];
        count1 =count1+1;
    end
    ha_count = [];
    ha_all1 = ha_all(:,1:8);
    for iii = 1:size(ha_all1,1)
        ha_count(iii) = sum(not(isnan(ha_all1(iii,:))));
    end
    ha_significance = sum(ha_all1,2,'omitnan')';
    sig_ratio{count2} = ha_significance./ha_count;
    flow_var_all{count2} = flow_var_max;
    pp_anova{count2} = pp_all;
    %%
    count2 =count2+1;
end
%%
count2  =1;
for current_area = [1,2,4]
    %%
    indx = find(contains(T.Area,area{current_area}));
    current_T = T(indx,:);
    %%
    list = 1:size(current_T,1); 
    flow_real = [];
    flow_permutation = [];
    count1 = 1;
    for kk = list
        ops = get_session_info2(current_T,kk,data_folder);
        fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
        ffname = fullfile(rfolder,fname);
        area_folder = area{current_area};
        load(fullfile(rfolder, fname));
        flow_var1 = 1-flow_var;
        %%
        flow_real(count1,:) = flow_var1(1,:);
        flow_permutation(count1,:) = mean(flow_var1(2:end,:),1);
        count1 = count1+1;
    end
    current_sig = sig_ratio{count2};
    sig_indx = (current_sig>=0);
    %%
    flow_real_all{count2} = flow_real(sig_indx,:);
    flow_permutation_all{count2} = flow_permutation(sig_indx,:);
    count2 = count2+1;
end
%%
color2 = {'g','r','c','m'};
h4i = figure('Renderer', 'painters', 'Position', [100 100 900 600]);
count1 = 1;
for current_area = [1,2,4]
    subplot(1,3,count1);
    flow_real_temp = flow_real_all{count1};
    flow_perm_temp = flow_permutation_all{count1};
    for k = 1:size(flow_real_temp,2)
        samplea = flow_real_temp(:,k);
        sample_n(k) = sum(not(isnan(samplea)));
    end
    flow_real_mean = mean(flow_real_temp,1,'omitnan');
    flow_real_std = std(flow_real_temp,1,'omitnan')./sqrt(sample_n);
    flow_perm_mean = mean(flow_perm_temp,1,'omitnan');
    flow_perm_std = std(flow_perm_temp,1,'omitnan')./sqrt(sample_n);
    errorbar(edges(1,1:8),flow_real_mean(1:8),flow_real_std(1:8),'color',color2{current_area});
    hold on;
    errorbar(edges(1,1:8),flow_perm_mean(1:8),flow_perm_std(1:8),'k');    
    xlim([0,0.01]);
    ylim([0,1]);
    yticks([0:0.2:1]);
    yticklabels(string(0:0.2:1));
    %% unbalanced two-way-anova
    flow_real_temp2 = flow_real_temp(:,1:8);
    flow_perm_temp2 = flow_perm_temp(:,1:8);
    flow_2way = [flow_real_temp2;flow_perm_temp2];
    edges2 = edges(1,2:size(flow_real_temp2,2)+1);
    edges_matrix2 = repmat(edges2,size(flow_2way ,1),1);
    permute_id_matrix2 = ones(size(flow_2way));
    permute_id_matrix2(1:size(flow_real_temp2,1),:) = 0;    
    flow_2way2 = flow_2way(:);    
    edges_matrix3 = edges_matrix2(:);
    permute_id_matrix3 = permute_id_matrix2(:);
    flow_2way3 = flow_2way2(not(isnan(flow_2way2)));
    edges_matrix3 = edges_matrix3(not(isnan(flow_2way2)));
    permute_id_matrix3  = permute_id_matrix3(not(isnan(flow_2way2)));
    pp2 = anovan(flow_2way3,{edges_matrix3 permute_id_matrix3},...
        'model',2,'varnames',{'amp','permute'},'display','off');
    pp_all_area(:,count1) = pp2;
    %%
    % indx1 = not(isnan(flow_var1(1,:)));
    h = []; p = [];
    for ii = 1:8
        flow_real_temp1 = flow_real_temp2(:,ii); 
        flow_perm_temp1 = flow_perm_temp2(:,ii); 
        [h(ii),p(ii)] = ttest(flow_real_temp1,flow_perm_temp1);
    end
    %%
    for ii  =1:8
        if not(isnan(h(ii)))
            if h(ii) 
                text(edges(1,ii),flow_real_mean(ii)+0.1,'*','fontsize',12);
            end
        end
    end
    count1 = count1+1;
end
%%
print(h4i, fullfile(save_folder,'Fig4i_flow_matchingIndex.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%%
count2  =1;
for current_area = [1,2,4]
    indx = find(contains(T.Area,area{current_area}));
    current_T = T(indx,:);
    list = 1:size(current_T,1);  
    phase_real = [];
    phase_permutation = [];
    count1 = 1;
    for kk = list
        ops = get_session_info2(current_T,kk,data_folder);
        fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
        ffname = fullfile(rfolder,fname);
        area_folder = area{current_area};
        load(fullfile(rfolder, fname));
        phase_var1 = 1-phase_var;
        %%
        phase_real(count1,:) = phase_var1(1,:);
        phase_permutation(count1,:) = mean(phase_var1(2:end,:),1);
        count1 = count1+1;
    end
    current_sig = sig_ratio{count2};
    sig_indx = (current_sig>=0);
    phase_real_all{count2} = phase_real(sig_indx,:);
    phase_permutation_all{count2} = phase_permutation(sig_indx,:);
    count2 = count2+1;
end
%%
color2 = {'g','r','c','m'};
h4h = figure('Renderer', 'painters', 'Position', [100 100 900 600]);
count1 = 1;
for current_area = [1,2,4]
    subplot(1,3,count1);
    phase_real_temp = phase_real_all{count1};
    phase_perm_temp = phase_permutation_all{count1};
    for k = 1:size(phase_real_temp,2)
        samplea = phase_real_temp(:,k);
        sample_n(k) = sum(not(isnan(samplea)));
    end
    phase_real_mean = mean(phase_real_temp,1,'omitnan');
    phase_real_std = std(phase_real_temp,1,'omitnan')./sqrt(sample_n);
    phase_perm_mean = mean(phase_perm_temp,1,'omitnan');
    phase_perm_std = std(phase_perm_temp,1,'omitnan')./sqrt(sample_n);
    errorbar(edges(1,1:8),phase_real_mean(1:8),...
        phase_real_std(1:8),'color',color2{current_area});
    hold on;
    errorbar(edges(1,1:8),phase_perm_mean(1:8),phase_perm_std(1:8),'k');    
    xlim([0,0.01]);
    ylim([0,1]);
    yticks([0:0.2:1]);
    yticklabels(string(0:0.2:1));
    %% unbalanced two-way-anova
    phase_real_temp2 = phase_real_temp(:,1:8);
    phase_perm_temp2 = phase_perm_temp(:,1:8);
    phase_2way = [phase_real_temp2;phase_perm_temp2];
    edges2 = edges(1,2:size(phase_real_temp2,2)+1);
    edges_matrix2 = repmat(edges2,size(phase_2way ,1),1);
    permute_id_matrix2 = ones(size(phase_2way));
    permute_id_matrix2(1:size(phase_real_temp2,1),:) = 0;    
    phase_2way2 = phase_2way(:);    
    edges_matrix3 = edges_matrix2(:);
    permute_id_matrix3 = permute_id_matrix2(:);
    phase_2way3 = phase_2way2(not(isnan(phase_2way2)));
    edges_matrix3 = edges_matrix3(not(isnan(phase_2way2)));
    permute_id_matrix3  = permute_id_matrix3(not(isnan(phase_2way2)));
    pp2 = anovan(phase_2way3,{edges_matrix3 permute_id_matrix3},...
        'model',2,'varnames',{'amp','permute'},'display','off');
    pp_all_area(:,count1) = pp2;
    
    h = []; p = [];
    for ii = 1:8
        phase_real_temp1 = phase_real_temp2(:,ii); 
        phase_perm_temp1 = phase_perm_temp2(:,ii); 
        [h(ii),p(ii)] = ttest(phase_real_temp1,phase_perm_temp1);
    end    
    for ii  =1:8
        if not(isnan(h(ii)))
            if h(ii) 
                text(edges(1,ii),phase_real_mean(ii)+0.1,'*','fontsize',12);
            end
        end
    end
    count1 = count1+1;
end
%%
print(h4h, fullfile(save_folder,'Fig4h_phase_matchingIndex.pdf'),...
    '-dpdf', '-bestfit', '-painters');