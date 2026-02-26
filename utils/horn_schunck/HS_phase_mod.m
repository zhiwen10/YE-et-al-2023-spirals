function [u, v] = HS_phase_mod(im1, im2, alpha, ite, uInitial, vInitial, displayFlow, displayImg)
% Horn-Schunck optical flow method 
% Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
% 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337
%
% Usage:
% [u, v] = HS(im1, im2, alpha, ite, uInitial, vInitial, displayFlow)
% For an example, run this file from the menu Debug->Run or press (F5)
%
% -im1,im2 : two subsequent frames or images.
% -alpha : a parameter that reflects the influence of the smoothness term.
% -ite : number of iterations.
% -uInitial, vInitial : initial values for the flow. If available, the
% flow would converge faster and hence would need less iterations ; default is zero. 
% -displayFlow : 1 for display, 0 for no display ; default is 1.
% -displayImg : specify the image on which the flow would appear ( use an
% empty matrix "[]" for no image. )
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008
% Rev: Jan 2009
alpha=1;
ite=100;
uInitial = zeros(size(im1(:,:,1)));
vInitial = zeros(size(im2(:,:,1)));
displayFlow=0;
displayImg=im1;
%% Convert images to grayscale
%%
% Set initial value for the flow vectors
u = uInitial;
v = vInitial;
% Estimate spatiotemporal derivatives
[fx, fy, ft] = computeDerivatives_mod(im1, im2);
% Averaging kernel
kernel_1=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
% Iterations
for i=1:ite
    % Compute local averages of the flow vectors
%     uAvg=angle(exp(i*conv2(u,kernel_1,'same')));
%     vAvg=angle(exp(i*conv2(v,kernel_1,'same')));
    uAvg = conv2(u,kernel_1,'same');
    vAvg = conv2(v,kernel_1,'same');
%     uAvg = conv2ang(u,kernel_1,'same');
%     vAvg = conv2ang(v,kernel_1,'same');
    % Compute flow vectors constrained by its local average and the optical flow constraints
    u= uAvg - ( fx .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2); 
    v= vAvg - ( fy .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
end

u(isnan(u))=0;
v(isnan(v))=0;

%% Plotting
if displayFlow==1
    plotFlow(u, v, displayImg, 5, 5); 
end