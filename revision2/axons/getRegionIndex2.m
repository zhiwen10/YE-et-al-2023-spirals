function st_region_indx = getRegionIndex2(region_label,st)
% region_label =  {"MOp", "MOs"}; % {"SSp"};
st_region_indx = [];
for r = 1:numel(region_label)
    st_region_indx = find(strcmp(st.acronym, region_label{r}));
end