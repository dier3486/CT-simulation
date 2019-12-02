function [detangle, eqangle] = detpos2angle(detposition, focalposition)
% [detangle, eqangle] = detpos2angle(detposition, focalposition)
% where the 'detposition' could be SYS.detector.position(1:Npixel, :);
% and the 'focalposition' could be SYS.source.focalposition(1, :);

Npixel = size(detposition, 1);
% angles
detangle = atan2(detposition(1:Npixel, 2) - focalposition(2), ...
    detposition(1:Npixel, 1) - focalposition(1));
midangle = atan2(-focalposition(2), -focalposition(1));

midindex =  find(detangle<midangle, 1, 'last');
midc = midindex + (midangle-detangle(midindex))/(detangle(midindex+1) - detangle(midindex));

% to equal angle
xx = (1:Npixel)' - midc;
delta = sum(xx.*(detangle-midangle))/sum(xx.^2);
eqangle = xx.*delta + midangle;
end