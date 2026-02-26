function explained_var_all = get_variance_explained(U,dV_raw,dV_predict)
Ur = double(reshape(U,size(U,1)*size(U,2),size(U,3)));
%% seperate prediction
for i = 1:size(dV_raw,3)
    traceE = Ur*dV_raw(:,:,i);
    traceEp = Ur*dV_predict(:,:,i);
    explained_var3 = sseExplainedCal(traceE ,traceEp);
    explained_var3 = reshape(explained_var3,size(U,1),size(U,2));
    explained_var_all(:,:,i) = explained_var3;
end