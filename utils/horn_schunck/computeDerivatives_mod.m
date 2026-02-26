function [fx, fy, ft] = computeDerivatives(im1, im2)

if size(im2,1)==0
    im2=zeros(size(im1));
end

% [gfx1, gfy1, badChannels1] = phasegradient(im1);
% [gfx2, gfy2, badChannels2] = phasegradient(im2);
% ft = anglesubtract(im1,im2);
% fx = gfx1+gfx2;
% fy = gfy1+gfy2;


% Horn-Schunck original method
% fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
% fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
% ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
% % 
% fx = angle(exp(i*fx));
% fy = angle(exp(i*fy));
% ft = angle(exp(i*ft));

% Horn-Schunck original method
% fx = angle(exp(i*conv2(im1,0.25* [-1 1; -1 1],'same'))) + angle(exp(i*conv2(im2, 0.25*[-1 1; -1 1],'same')));
% fy = angle(exp(i*conv2(im1, 0.25*[-1 -1; 1 1], 'same'))) + angle(exp(i*conv2(im2, 0.25*[-1 -1; 1 1], 'same')));
% ft = angle(exp(i*conv2(im1, 0.25*ones(2),'same'))) + angle(exp(i*conv2(im2, -0.25*ones(2),'same')));
% 
% fx = angle(exp(i*fx));
% fy = angle(exp(i*fy));
% ft = angle(exp(i*ft));
% derivatives as in Barron
% fx= conv2(im1,(1/12)*[-1 8 0 -8 1],'same');
% fy= conv2(im1,(1/12)*[-1 8 0 -8 1]','same');
% ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
% fx=-fx;fy=-fy;

% by Nick Steinmetz
% fx = angle(exp(1i*diff(im1,[],2)))/-2+ angle(exp(1i*diff(im2,[],2)))/-2;
% fy = angle(exp(1i*diff(im1,[],1)))/-2+ angle(exp(1i*diff(im2,[],1)))/-2;
% ft = angle(exp(1i*(im1-im2))); 

% by Nick Steinmetz
fx = angle(exp(i*conv2(im1, [-1 1], 'same')))/2+ angle(exp(i*conv2(im2, [-1 1], 'same')))/2;
fy = angle(exp(i*conv2(im1, [-1; 1], 'same')))/2+ angle(exp(i*conv2(im2, [-1; 1], 'same')))/2;
ft = angle(exp(i*(im1-im2))); 
% 
fx = angle(exp(i*fx));
fy = angle(exp(i*fy));

