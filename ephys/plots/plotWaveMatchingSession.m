function hs13f = plotWaveMatchingSession(T,data_folder,save_folder)
%%
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
hs13f = figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
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
            [h(ii),p(ii)] = ttest2(flow_var1(1,ii),flow_var1(2:end,ii),...
                "Tail","right");
        end
        ha = nan(1,20); ha(1,indx1) = h;
        pa = nan(1,20); pa(1,indx1) = p;
        pa_all = [pa_all;pa];
        ha_all = [ha_all;ha];
        flow_var_max = [flow_var_max; max(flow_var1(1,:))];       
        max_flow = max(flow_var1,[],1);
        %%
        pos = kk+(count2-1)*12;
        ax2 = subplottight(3,12,pos);
        ax2.Position(1) = ax2.Position(1)+0.018;
        ax2.Position(2) = ax2.Position(2)+0.03;
        ax2.Position(3) = ax2.Position(3)-0.025;
        ax2.Position(4) = ax2.Position(4)-0.08;
        scatter(edges1,flow_var1(1,:),8,color1{count2});
        hold on;
        plot(edges1,flow_var1(1,:),color1{count2});
        for i = 2:size(flow_var1,2)
            scatter(edges1,flow_var1(i,:),8,'k');
            hold on;
            plot(edges1,flow_var1(i,:),'k');
        end
        for i = 1:numel(edges1)
            if ha(i)==1 & edges1(i)<=0.01
                text(edges1(1,i)-0.0001,max_flow(i)+0.2,'*','fontsize',12);
            end
        end
        xlim([0,0.01]);
        ylim([0,1]);
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        xticks([0, 0.005, 0.01])
        xticklabels({'0','0.5','1'})
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
    
    count2 =count2+1;
end
%%
print(hs13f, fullfile(save_folder,'FigS13f_wave_MatchingIndex.pdf'),...
    '-dpdf', '-bestfit', '-painters');