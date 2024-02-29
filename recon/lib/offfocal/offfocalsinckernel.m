function kernel = offfocalsinckernel(offintensity, offwidth, offedge, offsample)
% subfunction of off-focal kernel, sinc style
% kernel = offfocalsinckernel(offintensity, offwidth, offedge, offsample)
% NOTE: the offwidth should be the scaled value: offwidth = corr.offwidth/SID/(max(t_off)-min(t_off));

% regule the shape
offintensity = offintensity(:)';
offwidth = offwidth(:)';
offedge = offedge(:)';

tt = [0:offsample/2, -offsample/2+1:-1]';
kernel = sinc(tt * offwidth).*offintensity;
% I know offwidth = offwidth/SID/(max(t_off)-min(t_off));

% edge smooth
edge_smooth = 1./(1+exp((-offedge+abs(tt)./(offsample/2)).*(10./(0.5-abs(offedge-0.5)))));
edge_smooth = fillmissing(edge_smooth, 'nearest');

kernel = sum(kernel.*edge_smooth, 2);

end