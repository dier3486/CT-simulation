function [Tz, Tchn, Weight] = helicalbpinterp(Zgrid_f, Zv, Eta, Zeta, Cd, Nprot, Nslice, Zscale, ...
    delta_d, midchannel, Nimage, ConeWeightScale, ViewRange, concyclic)
% the subfucntion used in helical-BP

% I know Zgrid_f = Zgrid - Zshift.

if nargin < 14
    concyclic = false;
    % concyclic-mode 
end

% Phi, Deta
Phi = asin(Eta)./(pi*2);
Deta0 = sqrt(1-Eta.^2);

% replace the Deta by 
if ~concyclic
    Detapi = Deta0 - Zeta;
    Deta0 = Deta0 + Zeta;
else
    Detapi = 1 - Zeta./Deta0;
    Deta0 = 1 + Zeta./Deta0;
end

% Tz
Tz = (Zgrid_f - (Zv/Nprot - Phi).*Cd)./Deta0;
Tz = Tz.*Zscale + (Nslice+1)/2;

% Cone weight
Wcone = ((Nslice-1)/2 - abs(Tz)).*ConeWeightScale;
Wcone = Wcone.*(Wcone>0) - (Wcone-1).*(Wcone>1);

% interp target on Eta
Tchn = repmat(Eta./delta_d + midchannel, 1, Nimage);

% % Zupsampling-times upsampling
% Tz = (Tz - 1).*Zupsampling + 1;

% I know
ConeWeightScale_Cz = ConeWeightScale*Cd*Zscale;
sigma_z = (Nslice-1)/2/Cd/Zscale;
ZgvCd = Zgrid_f./Cd - Zv/Nprot;

% to normalize the cone weight
% a0 & Wa0, to consider the first related view on same side
a0 = ZgvCd + Phi - Deta0.*sigma_z;
ceila0 = ceil(a0);
Wa0 = 1 - (ceila0 - a0)./Deta0.*ConeWeightScale_Cz;
% I know Wa0 shall >=0
Wa0 = Wa0.*(Wa0>0) + ceila0;
Ca0 = ceil((ViewRange(1)-Zv)/Nprot);  % Note: ViewRange, Zv and Nport are integers
% while (ceila0 < Ca0) Wa0 shall = Ca0
Wa0 = Wa0 + (Ca0-Wa0).*(ceila0 < Ca0);
% I know, only when the view is closing to start the ceila0 could <Ca0

% b0 & Wb0, to consider the last related view on same side
b0 = ZgvCd + Phi + Deta0.*sigma_z;
floorb0 = floor(b0);
Wb0 = 1 - (b0 - floorb0)./Deta0.*ConeWeightScale_Cz;
% I know Wb0 shall >=0
Wb0 = floorb0 - Wb0.*(Wb0>0);
Cb0 = floor((ViewRange(2)-Zv)/Nprot);
% while (floorb0 > Cb0) Wb0 shall = Cb0 + 1
Wb0 = Wb0 + (Cb0 - Wb0).*(floorb0 > Cb0) + 1;
% I know, only when the view is closing to the end the floorb0 could >Cb0.

% api & Wapi, to consider the first related view on opposite
api = ZgvCd - Phi - Detapi.*sigma_z - 1/2;
ceilapi = ceil(api);
Wapi = 1 - (ceilapi - api)./Detapi.*ConeWeightScale_Cz;
Wapi = Wapi.*(Wapi>0) + ceilapi;
Capi = ceil((ViewRange(1)-Zv)/Nprot - 1/2);
Wapi = Wapi + (Capi - Wapi).*(ceilapi < Capi);

% bpi & Wbpi, to consider the last related view on opposite
bpi = ZgvCd - Phi + Detapi.*sigma_z - 1/2;
floorbpi = floor(bpi);
Wbpi = 1 - (bpi - floorbpi)./Detapi.*ConeWeightScale_Cz;
Wbpi = floorbpi - Wbpi.*(Wbpi>0);
Cbpi = floor((ViewRange(2)-Zv)/Nprot - 1/2);
Wbpi = Wbpi + (Cbpi - Wbpi).*(floorbpi > Cbpi) + 1;

% normalization upon the above components
Weight = Wcone./(Wb0 - Wa0 + Wbpi - Wapi);

end