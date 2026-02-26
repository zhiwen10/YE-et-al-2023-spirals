function [unique_spirals1] = density_color_plot2(pwAllRaw,histbin)
%%
histbin_half = histbin/2;
unique_spirals1 = unique(pwAllRaw(:,1:2),'rows');
%%
for k = 1:size(unique_spirals1,1)
    clear indx
    current_spiral = unique_spirals1(k,:);
    indx = find(pwAllRaw(:,1)>=current_spiral(1,1)-histbin_half & pwAllRaw(:,1)<=current_spiral(1,1)+histbin_half ...
        & pwAllRaw(:,2)>=current_spiral(1,2)-histbin_half & pwAllRaw(:,2)<=current_spiral(1,2)+histbin_half);
    unique_spirals1(k,3) = numel(indx);
end
