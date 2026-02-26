function [vxRaw2,vyRaw2] = flow_vector_scale(vxRaw,vyRaw,skip,zoom_scale)
vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
% skip = 3;
% zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
end
