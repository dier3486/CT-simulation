function bowtiecurve = bowtiecurvesample(minthick, maxthick, length, shoulderwidth, shoulderslope)
% a sample to design a bowtie curve
% bowtiecurve = bowtiecurvesample(maxthick, minthick, length,
% shoulderwidth, shoulderslope);
% Don't known how to design a bowtie? try this

xx = 0:0.1:length/2;
xx = [-fliplr(xx(2:end)) xx];

x1 = log(exp(abs(xx)./shoulderwidth)-1).*shoulderwidth;
yy = (maxthick-minthick)./(1+exp((shoulderwidth-x1).*shoulderslope)) + minthick;

bowtiecurve = [xx(:) yy(:)];