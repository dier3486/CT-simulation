% ball z cali test

% load('E:\data\simulation\PG\ball1.mat');

Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
detector = prmflow.system.detector;
focalpos = prmflow.system.focalposition(1, :);
detector.position = reshape(detector.position, Npixel, Nslice, 3);
fanangle = atan2(detector.position(:,:,2)-focalpos(2), detector.position(:,:,1)-focalpos(1));
detz = detector.position(:,:,3);
SDD = double(detector.SDD);
SID = double(detector.SID);

dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
rawdata = dataflow.rawdata(:, 1:Nviewprot);

[maxval, maxindex] = max(rawdata, [], 1);

maxcut = 0.5;
maxslice = floor((maxindex-1)/Npixel);
maxpixel = mod(maxindex-1, Npixel) + 1;

s = (maxval>maxcut) & (maxslice>2) & (maxslice<Nslice-1) & (maxpixel>3) & (maxpixel<Npixel-2);

viewangle = double(dataflow.rawhead.viewangle(s));
viewindex = 1:Nviewprot;
viewindex = viewindex(s);

rawdata = reshape(rawdata, Npixel, Nslice*Nviewprot);
maxslice = maxslice(s);
vpindex = maxslice + (viewindex-1).*Nslice;
% vpindex = vpindex(s);
datapixel = rawdata(:, vpindex);

rawdata = reshape(permute(reshape(rawdata, Npixel, Nslice, Nviewprot), [2 1 3]), Nslice, Npixel*Nviewprot);
maxpixel = maxpixel(s);
vslindex = maxpixel + (viewindex-1).*Npixel;
% vslindex = vslindex(s);
dataslice = rawdata(:, vslindex);

cut = 0.03;
datapixel = datapixel - min(datapixel);
datapixel(datapixel<max(datapixel).*cut) = 0;
dataslice = dataslice - min(dataslice);
dataslice(dataslice<max(dataslice).*cut) = 0;


Nvs = sum(s);
cpixel = zeros(1, Nvs);
cslice = zeros(1, Nvs);
for iview = 1:Nvs
    % pixel
    y = datapixel(:, iview);
    [~, ipeak] = max(y);
    t1 = find(y(1:ipeak)==0, 1, 'last');
    if isempty(t1)
        t1 = 1;
    end
    t2 = find(y(ipeak:end)==0, 1 , 'first')+ipeak-1;
    if isempty(t2)
        t2 = Npixel;
    end
    y = [0; y(t1:t2); 0];
    x = fanangle(t1:t2,maxslice(iview));
    pcs = spline(x, y);
    cpixel(iview) = ppweightcenter(pcs);
    
    % slice
    y = [0; dataslice(:, iview); 0;];
    x = detz(maxpixel(iview), :);
    xx = linspace(x(1), x(end), Nslice*20);
    yy = spline(x, y, xx);
    ycut = max(y(2), y(end-1));
    yy(yy<ycut) = 0;
    cslice(iview) = (yy*xx(:))/sum(yy);
end

xp = SDD./tan(cpixel);
zp = cslice.*csc(cpixel);

v0 = [100, 0, 0, 0, 0 , 0];
v = lsqnonlin(@(x) ballfitfun(viewangle, SID, SDD, x, double([xp; zp])), v0);

pv = ballprojectfun_test1(viewangle, SID, SDD, v);
[Rb, Zb, phib, xdet , zdet, alphadet] = tac(v);

viewangle_b = mod(viewangle+phib, pi*2);
sclose =  (viewangle_b<pi/2) | (viewangle_b>pi*3/2);

p_err = zp(sclose) - pv(2, sclose);
x_err = xp(sclose);

pixel_err = interp1(x_err, p_err, SDD./tan(mean(fanangle,2))).*sin(mean(fanangle,2));
% almost

function r = ballfitfun(viewangle, SID, SDD, x, p0)

p = ballprojectfun_test1(viewangle, SID, SDD, x);

r = p(:) - p0(:);

end
