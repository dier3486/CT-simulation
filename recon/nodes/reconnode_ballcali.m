function [dataflow, prmflow, status] = reconnode_ballcali(dataflow, prmflow, status)
% 'ball' calibration
% [dataflow, prmflow, status] = reconnode_detectorcali(dataflow, prmflow, status)

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% parameters and data to use in prmflow
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

% parameters to from pipe
caliprm = prmflow.pipe.(status.nodename);


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
options = optimoptions('lsqnonlin', 'Display', 'off');
v = lsqnonlin(@(x) ballfitfun(viewangle, SID, SDD, x, double([xp; zp])), v0, [], [], options);

pv = ballprojectfun_test1(viewangle, SID, SDD, v);
[Rb, Zb, phib, xdet, zdet, alphadet] = tac(v);

viewangle_b = mod(viewangle+phib, pi*2);
sclose =  (viewangle_b<pi/2) | (viewangle_b>pi*3/2);

p_err =  - zp(sclose) + pv(2, sclose);
x_err = xp(sclose);

pixel_err = interp1(x_err, p_err, SDD./tan(mean(fanangle,2))).*sin(mean(fanangle,2));

% return
dataflow.ballcorr.focal_xerror = xdet;
dataflow.ballcorr.focal_zerror = zdet;
dataflow.ballcorr.detslope = -alphadet;
dataflow.ballcorr.pixel_err = pixel_err;

% disp
fprintf('\nBall calibration:\n');
fprintf('    %-20s %-20s %-20s\n', 'ball Z', 'focal error', 'slope');
fprintf('    %-20.2f %-20.2f %.2f¡ã\n\n', Zb, zdet, -alphadet*(180/pi));

% plot
figure;
subplot(2,1,1);
hold on
plot(xp, zp);
plot(pv(1,:),pv(2,:));
grid on;
h2  = subplot(2,1,2);
plot(pixel_err);
axis([1, Npixel, -0.1, 0.1]);
set(h2, 'XTick', 1:16:Npixel);
set(h2, 'XTickLabel', cellfun(@num2str, num2cell(1:Npixel/16), 'UniformOutput', false));
grid on;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function r = ballfitfun(viewangle, SID, SDD, x, p0)

p = ballprojectfun(viewangle, SID, SDD, x);

r = p(:) - p0(:);

end


function p = ballprojectfun(viewangle, SID, SDD, Vin)
% test function for ball porjection fitting

[Rb, Zb, phib, x0 , z0, alpha0] = tac(Vin);

x = Rb.*sin(viewangle + phib);
y = Rb.*cos(viewangle + phib);

x = x.*SDD./(SID+y);
z = (SDD*Zb)./(SID+y);

Aalpha = [cos(alpha0) -sin(alpha0); sin(alpha0) cos(alpha0)];
p = Aalpha*[x; z] + [x0; z0];

end
