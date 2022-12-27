function interptable = ZetaEta2TzTable(Nsample, maxFOV, reconD, SID, Nslice, gantrytilt, coneflag)
% Zeta Eta to Tz table used in 3D back projection
% interptable = ZetaEta2TzTable(Nsample, maxFOV, reconD, SID, Nslice, gantrytilt, coneflag);

if ~exist('gantrytilt', 'var')
    % while nargin < 6
    gantrytilt = 0;
end
if ~exist('coneflag', 'var')
    % while nargin < 7
    coneflag = 1;
end

% recon FOV
Ds = single(min(maxFOV, reconD)/SID);
% I know reconD = sqrt(sum((recon.FOV/2+abs(recon.center)).^2))*2;
% zeta, eta sampling
zeta_samp = linspace(-Ds/2, Ds/2, Nsample);
eta_samp = linspace(-Ds/2, Ds/2, Nsample);
% (zeta, eta) mesh
[zeta_s, eta_s] = meshgrid(zeta_samp, eta_samp);
% homeomorphic axial reconstruction geometry
switch coneflag
    case 1
        % cone beam, full shots looping
        t_samp1 = axialconehomeomorph(eta_s, zeta_s, Nslice, gantrytilt);
        Nfill0 = single(2);
        Nleft = -Nslice/2+1-Nfill0;
        t_samp1(t_samp1<Nleft) = Nleft;
        Nright = Nslice/2+Nfill0;
        t_samp1(t_samp1>Nright) = Nright;
        t_samp1 = t_samp1 + Nslice/2+Nfill0;
        % zz sampling
        zz_samp = single(-Nslice+1:Nslice);
    case 2
        % "chord" beam, "half" shots looping
        t_samp1 = axialchordhomeomorph(eta_s, zeta_s, Nslice);
        Nfill0 = single(0);
        t_samp1 = cat(3, t_samp1(:,:,1).*0, t_samp1(:,:, Nslice+1:end), t_samp1(:,:,1).*0+Nslice+1);
        t_samp1 = t_samp1 + Nfill0;
        % zz sampling
        zz_samp = single(0:Nslice+1);
    otherwise
        error('Illegal coneflag %d!', coneflag);
end

% to return
interptable.Nfill0 = Nfill0;
interptable.Zeta = zeta_samp;
interptable.Eta = eta_samp;
interptable.zz = zz_samp;
interptable.t = t_samp1;

end
