function [fanangles, focalangle] = detpos2fanangles(detposition, focalposition)
% to calculate the fan angles by detecor position and focal position
% [fanangles, focalangle] = detpos2fanangles(detposition, focalposition);
% This sub function will be called in rebin, geomotry calibrations and off-focal correction. 

% fan angles
y = detposition(:, 2) - focalposition(:, 2)';
x = detposition(:, 1) - focalposition(:, 1)';
fanangles = atan2(y, x);
% I know the size(fanangles) is [Npixel, Nfocal] to support DFS.

% focal angle(s)
focalangle = atan2(-focalposition(:, 2), -focalposition(:, 1))';

return