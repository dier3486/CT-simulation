function interptable = omiga4table(Coeffgamma, Nsample, maxFOV, reconD, SID, Nslice, gantrytilt, coneflag)
% prepare the table for Sobolev space linearized omiga-4-points interpoloation method
% interptable = omiga4table(x);

if nargin<7
    gantrytilt = 0;
end
if nargin<8
    coneflag = 1;
end
if length(Nsample)<2
    Nsample(2) = Nsample(1);
end

% recon FOV
Ds = single(min(maxFOV, reconD)/SID);
% I know reconD = sqrt(sum((recon.FOV/2+abs(recon.center)).^2))*2;
% zeta, eta sampling
zeta_samp = linspace(-Ds/2, Ds/2, Nsample(1));
eta_samp = linspace(-Ds/2, Ds/2, Nsample(1));
% (zeta, eta) mesh
[zeta_s, eta_s] = meshgrid(zeta_samp, eta_samp);
% homeomorphic axial reconstruction geometry
switch coneflag
    case 1
        % cone beam, full shots looping
        t_samp1 = axialconehomeomorph(eta_s, zeta_s, Nslice, gantrytilt);
        Nfill0 = single(4);
        Nleft = -Nslice/2+1-Nfill0+3/2;
        t_samp1(t_samp1<Nleft) = Nleft;
        Nright = Nslice/2+Nfill0-3/2;
        t_samp1(t_samp1>Nright) = Nright;
        t_samp1 = t_samp1 + Nslice/2+Nfill0;
        % zz sampling
        zz_samp = single(-Nslice+1:Nslice);
    case 2
        % "chord" beam, "half" shots looping
        t_samp1 = axialchordhomeomorph(eta_s, zeta_s, Nslice);
        Nfill0 = single(1);
        t_samp1 = cat(3, t_samp1(:,:,1).*0, t_samp1(:,:, Nslice+1:end), t_samp1(:,:,1).*0+Nslice+1);
        t_samp1 = t_samp1 + Nfill0 + 1;
        % zz sampling
        zz_samp = single(0:Nslice+1);
    otherwise
        error('Illegal coneflag %d!', coneflag);
end

% alpha-beta interp prepare
Nfourp = single(Nsample(2));
index_intp = 1:Nfourp+1;
alpha_intp = single(linspace(0, 1, Nfourp+1)');
beta_intp = 1/2-sqrt(1+alpha_intp.*(1-alpha_intp).*4)./2;
fourpoint = [(1+alpha_intp-beta_intp)./2  (alpha_intp+beta_intp)./2  ...
            (Coeffgamma(1)/4)./sqrt(1-alpha_intp.*(1-alpha_intp).*Coeffgamma(2))];

% to return
interptable.Nfill0 = Nfill0;
interptable.Zeta = zeta_samp;
interptable.Eta = eta_samp;
interptable.zz = zz_samp;
interptable.t = t_samp1;
interptable.Nfourp = Nfourp;
interptable.fourpointindex = index_intp;
interptable.fourpoint = fourpoint;
interptable.convL = single([-1 2 -1]);

end
