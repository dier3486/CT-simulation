function [dataflow, prmflow, status] = reconnode_watergoback(dataflow, prmflow, status)
% cali node, put in sausage and go back to pig
% [dataflow, prmflow, status] = reconnode_watergoback(dataflow, prmflow,
% status);
% use to find the ideal target for BH and nonlinear calibration.

% prm
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
% from rebin
Nreb = prmflow.rebin.Nreb;
mid_u = prmflow.rebin.midchannel;
delta_d = prmflow.rebin.delta_d;

% caliprm
caliprm = prmflow.pipe.(status.nodename);

% forward filter
if isfield(caliprm, 'filter')
    H1 = loadfilter(filter, Nreb, delta_d);
else
    % default filter
    H1 = filterdesign('hann', Nreb, delta_d, 1.0);
end
% backward filter
H0 = filterdesign('', Nreb, delta_d, inf);
Hlen = length(H1);

% fill up
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, Nview);
A = zeros(Hlen, Nslice, Nview);
for ii = 1:Nslice
    [A(:, ii, :), n_left] = translatefillup(squeeze(dataflow.rawdata(:, ii, :)), Hlen, mid_u);
end

% reshape
Nd = Nslice*Nview;
A = reshape(A, Hlen, Nd);

% filter by H1
A = ifft(fft(A).*H1, 'symmetric');

% target
C0 = 1000.*(2/pi);
% cut
s = A>C0.*0.9;
% set water's projection to C0
C1 = zeros(1, Nd);
ix = zeros(2, Nd);
m = 16;
for ii = 1:Nd
    ix(1, ii) = find(s(:,ii), 1, 'first') + m;
    ix(2, ii) = find(s(:,ii), 1, 'last') - m;
    C1(ii) = mean(A(ix(1,ii):ix(2,ii), ii), 1);
    A(ix(1,ii):ix(2,ii), ii) = C1(ii);
end
% inverse filter by H0
A = ifft(fft(A)./H0, 'symmetric');
% scale to C0
A = A.*(C0./C1);
% interation one more step
A = ifft(fft(A).*H1, 'symmetric');
for ii = 1:Nd
    A(ix(1,ii):ix(2,ii), ii) = C0;
end
% inverse filter by H0 again
A = ifft(fft(A)./H0, 'symmetric');

% to return
% dataflow.rawdata = A((1:Npixel)+n_left, :);
dataflow.rawdata = A;
% the returns for inverse rebin
prmflow.rebin.Nreb = Hlen;
prmflow.rebin.Nleft = n_left;
prmflow.rebin.midchannel = mid_u+n_left;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end