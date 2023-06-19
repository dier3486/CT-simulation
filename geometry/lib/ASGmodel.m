function [Nt, alphaL, Dasg] = ASGmodel(normV, ASGdata)
% ASG model
% [Nt, alphaL, Dasg] = ASGmodel(normV, ASGdata);
%   where normV = normr(V);  V is B-A , e.g. V = detposition - focalposition 
% 

nVdotndet = dot(normV, ASGdata.normvector, 2);

t = [ASGdata.nAdotnXondet - dot(normV, ASGdata.nVedgex, 2)./nVdotndet  ...
     ASGdata.nAdotnZondet - dot(normV, ASGdata.nVedgez, 2)./nVdotndet];
t = abs(t.*ASGdata.ASGheight - ASGdata.gap_n)./ASGdata.pixelspace;

Nt = floor(t);
alphaL = ((t-Nt).*ASGdata.pixelspace - ASGdata.gap_p)./(ASGdata.pixelspace - ASGdata.gap_p.*2);
alphaL(alphaL<0) = 0;
alphaL(alphaL>1) = 1;

Dasg = ASGdata.ASGthickness./abs([dot(normV, ASGdata.nASGplx, 2)  dot(normV, ASGdata.nASGplz, 2)]);

end