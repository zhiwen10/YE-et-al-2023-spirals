function [unique_spirals1,scolor,low_color_bound,high_color_bound] = density_color_plot(pwAllRaw,histbin,varargin)
%%
% varargin{1} = low_color_bound;
% varargin{2} = high_color_bound;
histbin_half = histbin/2;
%%
unique_spirals1 = unique(pwAllRaw(:,1:2),'rows');
%%
for k = 1:size(unique_spirals1,1)
    clear indx
    current_spiral = unique_spirals1(k,:);
    indx = find(pwAllRaw(:,1)>=current_spiral(1,1)-histbin_half & pwAllRaw(:,1)<=current_spiral(1,1)+histbin_half ...
        & pwAllRaw(:,2)>=current_spiral(1,2)-histbin_half & pwAllRaw(:,2)<=current_spiral(1,2)+histbin_half);
    unique_spirals1(k,3) = numel(indx);
end
%%    
if nargin ==3
    low_color_bound = varargin{1};
    high_color_bound = varargin{2};
else
    low_color_bound = min(unique_spirals1(:,3));
    high_color_bound = max(unique_spirals1(:,3));
end
rgbColor = remapColor(unique_spirals1(:,3),low_color_bound,high_color_bound);
rgbColor(rgbColor(:)<=0) = 1; % avoid index ==0
rgbColor(rgbColor(:)>=255) = 255; % avoid index ==0
% cmap1 = flipud(cbrewer('div', 'RdYlBu', 255,'linear'));
% cmap1(cmap1(:)>1) = 1; cmap1(cmap1(:)<0) = 0;
cmap1 = hot(255);
scolor = cmap1(rgbColor,:);
%%
% imgray = ones(size(projectedTemplate1));
% % imrgb = cat(3,imgray,imgray,imgray);
% figure;
% % im1 = imshow(imrgb,[]);
% im1 = imshow(imgray,[]);
% set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
% hold on;
% scatter(unique_spirals1(:,1)*2,unique_spirals1(:,2)*2,6,scolor,'filled');
% hold on;
% scale2=1;
% overlayOutlines(coords,scale2);
% set(gca,'Ydir','reverse')