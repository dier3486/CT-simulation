function Crs = offfocalzcrossmatrix(Nslice, crossrate)
% a subfunction of the off-focal correction to return the Z-cross matrix
% the % crossrate shall \in [0 1], which could be SDD*(SID-d)/SID/(SDD-d) where d is the distance from focal spot to 
% collimator blades, but configurable.
% Afix = Crs * A; 
% Afix shall in size (Nslice, []);


if crossrate==1
    Crs = eye(Nslice);
    return;
end

D = Nslice+1;
x = (1:Nslice) - (Nslice+1)/2;
x1 = x.*crossrate - D*(1-crossrate)/2;
x2 = x.*crossrate + D*(1-crossrate)/2;

z = (1:Nslice) - Nslice/2;
C1 = z - x1(:);
C2 = x2(:) - (z-1);
C1(C1>1) = 1;   C1(C1<0) = 0;
C2(C2>1) = 1;   C2(C2<0) = 0;

Crs = C1 + C2 - 1;
Crs = Crs./sum(Crs, 2);

end