a1 = normrnd(1,2,[1,100]);
b1 = normrnd(1.2,2,[1,100]);
values_all = [a1,b1];
diff_real  = mean(a1-b1);
count = numel(a1)+numel(b1);
count1 = numel(a1);
diff1 = [];
for perm = 1:10000
    tp1_indx = randperm(count, count1);
    tp2_indx = setdiff(1:count, tp1_indx);
    tp1_perm = values_all(tp1_indx);
    tp2_perm = values_all(tp2_indx);
    diff1(perm,1) = mean(tp1_perm-tp2_perm);    
end
%%
figure;
histogram(diff1);
hold on;
xline(diff_real,'r--')
%%
[h1,p1] = ttest(a1,b1);
%%
p_value_two_tailed = mean(abs(diff1) >= abs(diff_real));