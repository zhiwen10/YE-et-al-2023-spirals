%% get flow vector coordinates and color
[row, col] = find(not(isnan(vxRaw2)));
indx1 = find(not(isnan(vxRaw2)));
vx = vxRaw2(indx1);
vy = vyRaw2(indx1);

% get color for vectors from phasemap
% Your phase matrix (values from -pi to pi)
phaseData = framea; % your 2D matrix
% phaseData(~BW2) = nan;
% phaseData = reshape(phaseData, size(framea,1),size(framea,2));
% Normalize to [0, 1]
normData = (phaseData + pi) / (2 * pi);
% Convert to colormap indices (e.g., 256-level)
nColors = 256;
ind = round(normData * (nColors - 1)) + 1;
ind = min(max(ind, 1), nColors); % clamp
% Get the colormap (e.g., hsv, parula, jet, etc.)
cmap = colorcet('C06','N',nColors);
% % Convert to RGB image (H x W x 3)
rgbImage = ind2rgb(ind, cmap);
% rgbimage2 = imresize(rgbImage, [1320,1140]);
rgbImage_1d = reshape(rgbImage,[],3);
vector_color = rgbImage_1d(indx1,:);
figure;
for kk = 1:numel(row)
    imH1Raw4 = quiver(col(kk), row(kk),vx(kk),vy(kk),'color',vector_color(kk,:),'lineWidth',1,'autoScale','off');
    hold on;
end
set(gca,'YDir','reverse')
arrow_matrix = [col,row,vx,vy];
arrow_colors = vector_color;
save("flow_field_map.mat","arrow_matrix","arrow_colors");