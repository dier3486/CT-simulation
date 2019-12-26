% to test inverse filter
% I know d1 is rawdata after rebin
filtname = 'lann';
freqscale = 1.2;
Npixel = prmflow.recon.Npixel;
delta_d = prmflow.recon.delta_d;

H1 = filterdesign(filtname, Npixel, delta_d, freqscale);
len = length(H1);

a1 = d1(:,1);
a2 = a1;
a2(len) = 0;
a2 = ifft(fft(a2).*H1, 'symmetric');

C0 = 1000.*(2/pi);
ix1 = find(a2>C0,1,'first');
ix2 = find(a2>C0,1,'last');
m = 16;

C1 = mean(a2(ix1+m:ix2-m));
a3 = a2;
a3(ix1+m:ix2-m) = C1;

H0 = filterdesign('', Npixel, delta_d, 2);
b1 = ifft(fft(a3)./H0, 'symmetric');

b2 = b1.*(C0/C1);
b3 = ifft(fft(b2).*H1, 'symmetric');

b3(ix1+m:ix2-m) = C0;
b4 = ifft(fft(b3)./H0, 'symmetric');

b5 = ifft(fft(b4).*H1, 'symmetric');