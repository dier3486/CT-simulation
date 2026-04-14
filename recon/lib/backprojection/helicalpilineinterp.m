function [Tz, Tchn, Weight] = helicalpilineinterp(Zgrid_f, Zv, Eta, Zeta, Cd, Nprot, Nslice, Zscale, Zupsampling, ...
    delta_d, midchannel, Nimage, concyclic)
% pi-line BP

if nargin < 13
    concyclic = false;
    % concyclic-mode 
end

% Phi, Deta
Phi = asin(Eta)./(pi*2);
Deta = sqrt(1-Eta.^2);

if ~concyclic
    % Tz
    Tz = (Zgrid_f - (Zv/Nprot - Phi).*Cd)./(Deta+Zeta);
    % in Pi?
    PiC = (Tz.*Deta./Cd - Phi).*4;
else
    Tz = (Zgrid_f - (Zv/Nprot - Phi).*Cd).*Deta./(Deta+Zeta);
    PiC = (Tz./Cd - Phi).*4;
end
Weight = PiC>=-1 & PiC<1;

% Z scale
Tz = Tz.*Zscale + (Nslice+1)/2;
% extrap (for big pitch)
% Tz(Tz<1) = 1;  Tz(Tz>Nslice) = Nslice;
Tz = Tz.*(Tz>=1 & Tz<=Nslice) + (Tz<1).*1.0 + (Tz>Nslice).*Nslice;
% shift Tz to the Z-upsampled position
Tz = (Tz - 1).*Zupsampling + 1;

% interp target on Eta
Tchn = repmat(Eta./delta_d + midchannel, 1, Nimage);

end