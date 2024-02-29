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
Crs = offfocalzcrossmatrix(Nslice, crossrate);
% apply 
A = Crs*A;
% permute back
A = permute(reshape(A, Nslice, Npixel, Nview), [2 1 3]);
% done

end