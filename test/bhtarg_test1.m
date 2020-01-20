load('E:\data\rawdata\bhtest\A1_300.mat');

[Hlen, Nview] = size(A1);
% A2 = ifft(fft(A1).*H1, 'symmetric');

xx = 1:Hlen;
Cw = (xx*A1)./sum(A1);
m = 30;
Cw = smooth([Cw(end-m+1:end) Cw Cw(1:m)], 0.03, 'rloess');
Cw = Cw(m+1:end-m);

A2 = ifft(fft(A1).*H1, 'symmetric');

Cmid = mean(Cw);
Am = zeros(size(A1));
for ii = 1:Nview
    Am(:, ii) = interp1(xx, A2(:, ii), xx+(Cw(ii)-Cmid), 'linear', 0);
end
Ac = mean(Am,2);

C0 = 1000.*(2/pi);
x1_l = find(Ac>C0, 1, 'first');
x1_r = find(Ac>C0, 1, 'last');
m = 16;

Cmn = mean(Ac(x1_l+m:x1_r-m));
x2_l = find(Ac(x1_l+2:end)<Cmn, 1, 'first')+x1_l+1;
x2_r = find(Ac(1:x1_r-2)<Cmn, 1, 'last');

span = 50.0/(x2_r-x2_l);
Asm = smooth(Ac(x2_l:x2_r), span, 'rloess');
x3_l = find(diff(Asm>Cmn)<0, 1, 'first')+1;
x3_r = find(diff(Asm>Cmn)>0, 1, 'last');
Asm(x3_l:x3_r) = Cmn;

ix = zeros(2, Nview);
ix(1,:) = round(Cw-Cmid+x2_l);
ix(2,:) = round(Cw-Cmid+x2_r);

A3 = A2;
C1 = zeros(1, Nview);
for ii = 1:Nview
    xx = (x2_l:x2_r)+(Cw(ii)-Cmid);
    C1(ii) = mean(A3(ix(1,ii):ix(2,ii), ii));
    tofill = interp1(xx, Asm, ix(1,ii):ix(2,ii), 'linear', 'extrap');
    A3(ix(1,ii):ix(2,ii), ii) = tofill.*(C1(ii)/mean(tofill));
    C1(ii) = C1(ii)./mean(tofill)*Cmn;
end

% inverse filter by H0
A4 = ifft(fft(A3)./H0, 'symmetric');
% scale to C0
A4 = A4.*(C0./C1);
% interation one more step
A5 = ifft(fft(A4).*H1, 'symmetric');
% inverse filter by H0 again
A6 = ifft(fft(A5)./H0, 'symmetric');
