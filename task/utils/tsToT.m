function t = tsToT(ts, numSamps)
% function t = tsToT(ts, numSamps)
% 
% Convert an alf 'timestamps' representation into explicit times of each
% sample

t = interp1(ts(:,1), ts(:,2), [0:(numSamps-1)]');