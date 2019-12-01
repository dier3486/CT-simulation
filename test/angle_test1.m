% I know the SYS has been load

Npixel = double(SYS.detector.Npixel);
focalposition = double(SYS.source.focalposition(1,:));
detposition = double(SYS.detector.position);
% angles
detangles = atan2(detposition(1:Npixel, 2) - focalposition(2), ...
    detposition(1:Npixel, 1) - focalposition(1));
midangle = atan2(-focalposition(2), -focalposition(1));

midindex =  find(detangles<midangle, 1, 'last');
midc = midindex + (midangle-detangles(midindex))/(detangles(midindex+1) - detangles(midindex));
xx = (1:Npixel)' - midc;
% to equal angle
delta1 = mean(diff(detangles));
eqangle1 = xx.*delta1 + midangle;

delta2 = sum(xx.*(detangles-midangle))/sum(xx.^2);
eqangle2 = xx.*delta2 + midangle;