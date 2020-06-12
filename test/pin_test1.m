% pin cali main, test

%% #1 load data and air corr
pinxml = 'F:\data-Dier.Z\PG\bay3_0601\Pinrecon.xml';
pindata = 'F:\data-Dier.Z\PG\bay3_0601\pin\pin_1';
[~, dataflow_all, prmflow_all] = CRISrecon(pinxml, pindata);

%% #2 get 1st rot
Nview_all = prmflow_all.protocol.viewnumber;
Nviewprot = prmflow_all.protocol.viewperrot;
Nrot_all = Nview_all/Nviewprot;

Nrot = Nrot_all;
viewstep = Nviewprot/4;
Nstep = 20;

istep = 1;

startview = (istep-1)*viewstep + 1;
endview = startview + Nviewprot*Nrot - 1;

dataflow = dataflow_all;
dataflow.rawdata = dataflow.rawdata(:, startview:endview);
for ifield = fieldnames(dataflow.rawhead)'
    dataflow.rawhead.(ifield{1}) = dataflow.rawhead.(ifield{1})(:, startview:endview);
end

prmflow = prmflow_all;
prmflow.recon.Nview = Nviewprot*Nrot;

%% #3 prms and spline fitting
% inputs are dataflow and prmflow
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nrot = Nview/Nviewprot;
detector = prmflow.system.detector;
focalposition = prmflow.system.focalposition(prmflow.recon.focalspot, :);
Npixelpermod = 16;
Nmod = Npixel/Npixelpermod;
SID = double(detector.SID);
SDD = double(detector.SDD);
hz = double(detector.hz_ISO);
viewangle = double(dataflow.rawhead.viewangle);

% det to fanangle
[fanangles0, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles0 = fanangles0 - focalangle;
fanangles0 = reshape(fanangles0, Npixel, Nslice);

% ini
fanangles = fanangles0;

% splines
cs = pinsplines(dataflow.rawdata, fanangles, Npixel, Nslice, Nview);

%% #4 pin-model fitting
p0 = reshape(double([cs(:).p]), Nslice, Nview);

x0 = [200/SID  0  0       0       0         0     0       0       0       0];
%    [r,       phi,  zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, zshift]

% x0 = pf1;

% pinfit = lsqnonlin(@(x) pinfitfun(viewangle, Nslice-2, hz, x, p0(2:end-1, :)), x0);
pinfit = lsqnonlin(@(x) pinfitfun(viewangle, Nslice, hz, x, p0), x0);
p1 = pinfitfun(viewangle, Nslice, hz, pinfit, 0);
p1 = p1(1:Nslice, 1:Nview);
dp = p1 - p0;

%% #5 pin matrix
alpha_L = [1.0, 0.5];
[A, L, indexrange] = pinleftoptmatrix(cs, Npixel, Nslice, Nview, Npixelpermod, alpha_L);



%% pin fit subfunction
function r = pinfitfun(viewangle, Nslice, hz, x, p0)

if ~iscell(x)
    x = num2cell(x);
end

p = pinprojectfun(viewangle, Nslice, hz, x{:});

r = p - p0;
er = (r - mean(r)).*(Nslice-1);
% er = std(r).*Nslice;
dr = r(:,end) - r(:, 1);
r = [r er dr.*size(p, 2)];
% r = [r dr.*size(p, 2);
%      er 0];

end


