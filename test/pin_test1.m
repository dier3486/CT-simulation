
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nview = prmflow.recon.Nview;
detector = prmflow.system.detector;
focalposition = prmflow.system.focalposition(prmflow.recon.focalspot, :);

[fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles = fanangles - focalangle;
fanangles = reshape(fanangles, Npixel, Nslice);

dataflow.rawdata = reshape(dataflow.rawdata, Npixel, []);

[pmax, idxmax] = max(dataflow.rawdata, [], 1);

m = 20;
index = (-m:m)' + idxmax;
sout = (index>Npixel) | (index<1);
index(index>Npixel) = Npixel;
index(index<1) = 1;
index = index + (0:Nview*Nslice-1).*Npixel;
A1 = dataflow.rawdata(index);
A1(sout) = 0;
C1 = (sum(dataflow.rawdata, 1) - sum(A1, 1))./(Npixel-m*2-1-sum(sout, 1));

A1 = A1 - C1;
A1(sout) = 0;
pmax = pmax - C1;

cut = 0.05;
A1(A1 < repmat(pmax.*cut, m*2+1, 1)) = 0;

% ini cs
cs = spline([0 1], [0 0]);
cs.index = [];
cs.p = nan;
cs.dp = nan;
cs(Nview*Nslice) = cs(1);
for ii = 1:Nview*Nslice
    islice = mod(ii-1, Nslice) + 1;
    find_ii = find(A1(:, ii)>0, 1, 'first')-1 : find(A1(:, ii)>0, 1, 'last')+1;
    fanindex_ii = mod(index(find_ii, ii)-1, Npixel)+1;
    fan_ii = fanangles(fanindex_ii, islice);
    tmp = spline(fan_ii, [0; A1(find_ii, ii); 0]);
    tmp.index = fanindex_ii;
    [tmp.dp, tmp.p] = ppcenterderivative(tmp);
    cs(ii) = tmp;
end

