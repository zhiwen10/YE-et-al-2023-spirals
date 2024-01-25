% return unbalanced 2-way anova p values [pp]
% and p value for each bin with student t-test
%%
function [pp,p] = anova_test2(edges1,ratio_all,ratio_all_perm)
%% unbalanced two-way-anova
edges_matrix = repmat(edges1,size(ratio_all,1),1);
edges_matrix_all = [edges_matrix; edges_matrix];

id_matrix = ones(size(ratio_all));
id_matrix_permute = zeros(size(ratio_all_perm));
id_matrix_all = [id_matrix; id_matrix_permute];

ratio_matrix_all = [ratio_all; ratio_all_perm];
%% linearize 
ratio_matrix_all2 = ratio_matrix_all(:);
edges_matrix_all2 = edges_matrix_all(:);
id_matrix_all2 = id_matrix_all(:);
% indx = not(isnan(ratio_matrix_all2));
% ratio_matrix_all2 = ratio_matrix_all2(indx);
% edges_matrix_all2 = edges_matrix_all2(indx);
% id_matrix_all2 = id_matrix_all2(indx);
pp = anovan(ratio_matrix_all2,{edges_matrix_all2 id_matrix_all2},'model',2,'varnames',{'amp','permute'},'display','off');
%%
for ii = 1:size(ratio_all,2)
    [h(ii),p(ii)] = ttest(ratio_all(:,ii),ratio_all_perm(:,ii));
end