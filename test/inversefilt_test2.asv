% to test inverse filter
% I know d1 is rawdata after rebin
filtname = 'hann';
freqscale = 1.2;
Npixel = prmflow.recon.Npixel;
mid_u = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d;
Nview = prmflow.recon.Nview;

% forward filter
H1 = filterdesign(filtname, Npixel, delta_d, freqscale);
% backward filter
H0 = filterdesign('', Npixel, delta_d, inf);

len = length(H1);
n_l = floor((len+2-Npixel)/2);
n_r = ceil((len+2-Npixel)/2);

d1=dataflow.rawdata(:, 30:64:end);

% find model
xx = 1:Npixel;
[~, s_mid] = min(abs(xx*d1./sum(d1)-mid_u));
a0 = d1(:,s_mid);

% fillup
a1 = zeros(len, Nview);
a1(n_l:len-n_r+1, :) = d1;

Nfill = 100;
af_l = zeros(Nfill, Nview);
af_r = zeros(Nfill, Nview);
for ii = 1:Nview
    % R
    index_r = len-n_r;
    x_r = find(a1(:, ii)>a1(index_r, ii), 1, 'first');
    if x_r>n_r+1
        xx = x_r-1+(a1(index_r, ii)-a1(x_r-1, ii))/(a1(x_r, ii)-a1(x_r-1, ii))-(1:Nfill);
        af_r(:, ii) = interp1(a1(:, ii), xx);
    end
    % L
    index_l = n_l+1;
    x_l = find(a1(:, ii)>a1(index_l, ii), 1, 'last');
    if x_l<len-n_r
        xx = x_l+(a1(index_l, ii)-a1(x_l, ii))/(a1(x_l+1, ii)-a1(x_l, ii))+(Nfill:-1:1);
        af_l(:, ii) = interp1(a1(:, ii), xx);
    end
end
a1(n_l-Nfill+1:n_l, :) = af_l;
a1(len-n_r+1:len-n_r+Nfill, :) = af_r;

% A = a1;
% Nd = Nview;
% 
% % filter by H1
% a2 = ifft(fft(a1).*H1, 'symmetric');
% % target
% C0 = 1000.*(2/pi);
% % cut
% s = a2>C0.*0.9;
% % set water's projection to C0
% C1 = zeros(1, Nd);
% ix = zeros(2, Nd);
% m = 16;
% for ii = 1:Nd
%     ix(1, ii) = find(s(:,ii), 1, 'first') + m;
%     ix(2, ii) = find(s(:,ii), 1, 'last') - m;
%     C1(ii) = mean(a2(ix(1,ii):ix(2,ii), ii), 1);
%     a2(ix(1,ii):ix(2,ii), ii) = C1(ii);
% end
% % inverse filter by H0
% a3 = ifft(fft(a2)./H0, 'symmetric');
% % scale to C0
% a3 = a3.*(C0./C1);
% % interation one more step
% a4 = ifft(fft(a3).*H1, 'symmetric');
% for ii = 1:Nd
%     a4(ix(1,ii):ix(2,ii), ii) = C0;
% end
% % inverse filter by H0 again
% a5 = ifft(fft(a4)./H0, 'symmetric');
