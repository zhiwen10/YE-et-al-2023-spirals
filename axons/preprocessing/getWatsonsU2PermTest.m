function [p] = getWatsonsU2PermTest(T1)
[angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogram(T1);
for i = 1:100
    [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end