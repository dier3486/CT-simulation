function A = offfocalzcross(A, crossrate)
% a subfunction of the off-focal correction
% A = offfocalzcross(A, zrate); 
% A shall in shape (Npixel, Nslice, Nview);
% zrate \in [0 1], which should be SDD*(SID-d)/SID/(SDD-d) where d is the distance from focal spot to collimator blades

% size
[Npixel, Nslice, Nview] = size(A);
% permute
A = reshape(permute(A, [2 1 3]), Nslice, []);
% cross matrix
Crs = crossmatrix(Nslice, crossrate);
% apply 
A = Crs*A;
% permute back
A = permute(reshape(A, Nslice, Npixel, Nview), [2 1 3]);
% done

end

function Crs = crossmatrix(Nslice, zrate)

if zrate==1
    Crs = eye(Nslice);
    return;
end

D = Nslice+1;
x = (1:Nslice) - (Nslice+1)/2;
x1 = x.*zrate - D*(1-zrate)/2;
x2 = x.*zrate + D*(1-zrate)/2;

z = (1:Nslice) - Nslice/2;
C1 = z - x1(:);
C2 = x2(:) - (z-1);
C1(C1>1) = 1;   C1(C1<0) = 0;
C2(C2>1) = 1;   C2(C2<0) = 0;

Crs = C1 + C2 - 1;
Crs = Crs./sum(Crs, 2);

end