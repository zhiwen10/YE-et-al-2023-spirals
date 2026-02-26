function [Axy_all] = compare_flowfield_angle(data_raw,data_predict)
% this function calculate the angle difference between 2 flow vectors
% result in complex value
% function returns new vector, with angle representing angle difference (match is
% zero degrees),
% vector length represent raw flow vector length (higher value with higher
% travelling wave)
%%
x_size = size(data_raw.vxRaw, 2);
y_size = size(data_raw.vxRaw, 3);
frameN = size(data_raw.vxRaw, 1);
%%
Axy_all = [];
for i = 1:frameN    
    vxRaw =squeeze(data_raw.vxRaw(i,:,:)); 
    vyRaw = squeeze(data_raw.vyRaw(i,:,:));
    % vxRaw(BW) =  nan; vyRaw(BW) =  nan; 
    % [vxRaw1,vyRaw1] = flow_resize(vxRaw,vyRaw,skip);
    %%
    vxPredict =squeeze(data_predict.vxRaw(i,:,:)); 
    vyPredict = squeeze(data_predict.vyRaw(i,:,:));
    % vxPredict(BW) =  nan; vyPredict(BW) =  nan; 
    % [vxPredict1,vyPredict1] = flow_resize(vxPredict,vyPredict,skip);
    %% vectorize flow field in 2 cloumns
    vector1 = [vxRaw(:),vyRaw(:)];
    vector2 = [vxPredict(:),vyPredict(:)];
    %% get vector angle for each pixel
    A = vectorAngle(vector1, vector2);
    % A(A>pi) = A(A>pi)-2*pi;
    A2d = reshape(A,[x_size,y_size]);
    %% construct new vector 
    % cos(a) + i*sin(a) = exp(i*a)
    Axy = exp(1i*A2d);
    Axy_all(i,:,:) = Axy;
end