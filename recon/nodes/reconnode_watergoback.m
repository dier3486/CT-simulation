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

% smooth & fill
if isfield(caliprm, 'span')
    smooth_span = caliprm.span;
else
    smooth_span = 30;
end
% off-focal tolerance
if isfield(caliprm, 'offfocal')
    offfocal_tol = caliprm.offfocal;
else
    offfocal_tol = 0;
end
% to plot 
if isfield(caliprm, 'offplot')
    offplot = caliprm.offplot;
else
    offplot = false;
end

% forward filter
if isfield(caliprm, 'filter')
    H1 = loadfilter(caliprm.filter, Nreb, delta_d);
else
    % default filter
    H1 = filterdesign('hann', Nreb, delta_d, 1.0);
end
Hlen = length(H1);

% backward filter
H0 = filterdesign('', Nreb, delta_d, inf);

% fill up
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, Nview);
A = zeros(Hlen, Nslice, Nview);
if isfield(dataflow.rawhead, 'refblock')
    refblock = dataflow.rawhead.refblock;
else
    refblock = false(2, Nview);
end
for ii = 1:Nslice
    [A(:, ii, :), n_left] = translatefillup(squeeze(dataflow.rawdata(:, ii, :)), Hlen, refblock);
end

% reshape
Nd = Nslice*Nview;
A = reshape(A, Hlen, Nd);

% weight center
xx = 1:Hlen;
Wcnt = (xx*A)./sum(A, 1);
% smooth
Wcnt = reshape(Wcnt, Nslice, Nview);
m = 30;
for islice = 1:Nslice
    tmp = smooth([Wcnt(islice, end-m+1:end) Wcnt(islice, :) Wcnt(islice, 1:m)], 0.03, 'rloess');
    Wcnt(islice, :) = tmp(m+1:end-m);
end
Wmid = mean(Wcnt, 2);

% filter by H1
A = ifft(fft(A).*H1, 'symmetric');

% target
C0 = 1000.*(2/pi);

% mean shape & cut
A = reshape(A, Hlen, Nslice, Nview);
Amean = zeros(Hlen, Nslice);
xrange_l = zeros(Nslice, 1);
xrange_r = zeros(Nslice, 1);
Cmean = zeros(1, Nslice);
% index of not blocked
index_unblk = find(~any(refblock));
Nview_unblk = size(index_unblk(:),1);
for islice = 1:Nslice
    Ac = zeros(Hlen, Nview_unblk);
    xx = 1:Hlen;
    for iview_ub = 1:Nview_unblk
        iview = index_unblk(iview_ub);
        Ac(:, iview_ub) = interp1(xx, A(:, islice, iview), xx+(Wcnt(islice, iview)-Wmid(islice)), 'linear', 0);
    end
    Amean(:, islice) = mean(Ac, 2);
    % cut for mean
    x1_l = find(Amean(:, islice)>C0.*0.9, 1, 'first');
    x1_r = find(Amean(:, islice)>C0.*0.9, 1, 'last');
    m = 16;
    Cmean(islice) = mean(Amean(x1_l+m:x1_r-m, islice));
end

% to find water edge
Amean = mean(Amean, 2);
DfA = diff(Amean);
% left
[~, edge_pl] = max(Amean(x1_l: x1_l+m));
edge_pl = max(edge_pl, find(Amean(x1_l: x1_l+m)>C0.*1.3, 1, 'last'));
edge_pl = x1_l + edge_pl - 1;
xcut_l = find(DfA(edge_pl:end)>0, 1, 'first') + edge_pl - 1;
% right
[~, edge_pr] = max(Amean(x1_r-m: x1_r));
edge_pr = min(edge_pr, find(Amean(x1_r-m: x1_r)>C0.*1.3, 1, 'first'));
edge_pr = x1_r-m + edge_pr - 1;
xcut_r = find(DfA(1:edge_pr)<0, 1, 'last')+1;
%     % with off-focal
%     xcut_l(ii) = find(Amean(x1_l+2:end, ii)<Cmean(ii), 1, 'first')+x1_l+1;
%     xcut_r(ii) = find(Amean(1:x1_r-2, ii)<Cmean(ii), 1, 'last');    

A = reshape(A, Hlen, Nd);

% Asmth = zeros(xcut_r-xcut_l+1, 1);
span_smth = smooth_span/(xcut_r-xcut_l);
Asmth = smooth(Amean(xcut_l:xcut_r), span_smth, 'loess');
switch offfocal_tol
    case {'deep', 0}
        % deep off-focal
        x3_l = find(diff(Asmth>mean(Cmean))<0, 1, 'first')+1;
        x3_r = find(diff(Asmth>mean(Cmean))>0, 1, 'last');
        % with draw the cuts
        d_l = find(Asmth>mean(Cmean), 1, 'first');
        xrange_l(:) = xcut_l+d_l-1;
        d_r = find(Asmth>mean(Cmean), 1, 'last');
        xrange_r(:) = xcut_l+d_r-1;
        % fill Asmth
        Afill = repmat(Asmth, 1, Nslice);
        for ii = 1:Nslice
            Afill(:, ii) = Afill(:, ii).*(Cmean(ii)/mean(Cmean));
            Afill(x3_l:x3_r, ii) = Cmean(ii);
        end
    case {'weak', 1}
        % week off-focal
        x3_l = find(diff(Asmth>mean(Cmean))~=0, 1, 'first')+1;
        [~, x3_l] = max(Asmth(1:x3_l));
        x3_r = find(diff(Asmth>mean(Cmean))~=0, 1, 'last');
        [~, tmp] = max(Asmth(x3_r:end));
        x3_r = x3_r+tmp-1;
        % with draw the cuts
        xrange_l(:) = xcut_l+x3_l-1;
        xrange_r(:) = xcut_l+x3_r-1;
        % fill Asmth
        Afill = repmat(Cmean, xcut_r-xcut_l+1, 1);
    case {'none', 3}
        % no off-focal (sure?)
%         x3_l = 1;
%         x3_r = xcut_r-xcut_l+1;
        xrange_l(:) = xcut_l;
        xrange_r(:) = xcut_r;
        % replace Asmth
        Afill = repmat(Cmean, xcut_r-xcut_l+1, 1);
    otherwise
        error('Illegal off-focal tolerance: %s', num2str(offfocal_tol));
end

if offplot
    fg1 = figure;
    set(fg1, 'Position', [400 458 1000 420]);
    hold on;
    plot(Amean, 'b');
    a_plot = nan(size(Amean));
    a_plot(xcut_l:xcut_r) = Asmth;
    plot(a_plot, 'r');
    plot([xrange_l(1) xrange_r(1)], a_plot([xrange_l(1) xrange_r(1)]), 'r*');
    a_plot(xcut_l:xcut_r) = mean(Afill, 2);
    plot(a_plot, 'g');
    [tt, yt] = watersmile(Amean, C0, delta_d);
    plot(tt, yt);
    axis([xcut_l-20 xcut_r+20 mean(Cmean)-50 mean(Cmean)+50]);
    grid on;
    drawnow;
end

% index range of the target
index_range= zeros(2, Nd);
index_range(1, :) = reshape(round(Wcnt-Wmid+xrange_l), 1, []);
index_range(2, :) = reshape(round(Wcnt-Wmid+xrange_r), 1, []);

% replace
C1 = zeros(1, Nd);
for ii = 1:Nd
    % index
    islice = mod(ii-1, Nslice)+1;
    ix1 = index_range(1,ii);  
    ix2 = index_range(2,ii);
    % interp
    xx = (xcut_l : xcut_r) + Wcnt(ii)-Wmid(islice);
    index_sm = 1:xcut_r-xcut_l+1;
    tofill = interp1(xx, Afill(index_sm, islice), ix1:ix2, 'linear', 'extrap');
    % norm and replace
    C1(ii) = mean(A(ix1:ix2, ii))/mean(tofill);
    A(ix1:ix2, ii) = tofill.*C1(ii);
    C1(ii) = C1(ii)*Cmean(islice);
end

% inverse filter by H0
A = ifft(fft(A)./H0, 'symmetric');
% scale to C0
A = A.*(C0./C1);

% % interation one more step (skip)
% A = ifft(fft(A).*H1, 'symmetric');
% for ii = 1:Nd
%     A(index_range(1,ii):ix(2,ii), ii) = C0;
% end
% % inverse filter by H0 again
% A = ifft(fft(A)./H0, 'symmetric');

% % regular by max (skip)
% A = reshape(A, Hlen, Nslice, Nview);
% for ii = 1:Nslice
%     maxA = max(A(:,ii,:),[],1);
%     A(:,ii,:) = A(:,ii,:)+(mean(maxA)-maxA);
% end

% reshape
A = reshape(A, Hlen*Nslice, Nview);
index_range = reshape(index_range, 2*Nslice, Nview);

% to return
dataflow.rawdata = A;
dataflow.rawhead.index_range = index_range;
% the returns for inverse rebin
% prmflow.rebin.Nreb = Hlen;
% prmflow.rebin.Nleft = n_left;     % deleted
prmflow.recon.Npixel = Hlen;
prmflow.recon.midchannel = mid_u+n_left;
% NOTE: do no use prmflow.rebin to return values!

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [tt, yt] = watersmile(Amean, C0, delta_d)

tshell = 5;
tscale = 1.2;
delta_d = double(tshell/delta_d*(tscale-1));

n1 = find(Amean>C0, 1, 'first');
n2 = find(Amean>C0, 1, 'last');

d1 = n1 - (Amean(n1) - C0)/(Amean(n1) - Amean(n1-1)) - delta_d;
d2 = n2 + (Amean(n2) - C0)/(Amean(n2) - Amean(n2+1)) + delta_d;

d0 = (d1+d2)/2;
r0 = (d2-d1)/2;

m = 16;
xx = find(Amean(n1+m : n2-m) < C0) + n1 + m - 1;
tt = [d1 n1:n2 d2];

options = optimoptions('lsqnonlin','Display','off');
a = lsqnonlin(@(a) (1+(xx-d0)./r0).*(1-(xx-d0)./r0).*(a(1)+a(2).*((xx-d0)./r0).^2) - (Amean(xx)-C0), [0.1 0], [], [], options);
yt = C0 + (1+(tt-d0)./r0).*(1-(tt-d0)./r0).*(a(1)+a(2).*((tt-d0)./r0).^2);

end