function st_region_indx = getRegionIndex(region_label,st)
% region_label =  {"MOp", "MOs"}; % {"SSp"};
st_region_indx = [];
for r = 1:numel(region_label)
    region_id = st.id(strcmp(st.acronym, region_label{r}));
    st_region_indx = [st_region_indx; find(cellfun(@(x)contains(x, sprintf('/%d/', region_id)), st.structure_id_path))];
end