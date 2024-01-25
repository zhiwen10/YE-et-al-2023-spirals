function explained_var = get_variance_explained2(U,dV_raw,dV_predict)
Ur = double(reshape(U,size(U,1)*size(U,2),size(U,3)));
%% seperate prediction
traceE = Ur*dV_raw;
traceEp = Ur*dV_predict;
explained_var = sseExplainedCal(traceE ,traceEp);
explained_var = reshape(explained_var,size(U,1),size(U,2));