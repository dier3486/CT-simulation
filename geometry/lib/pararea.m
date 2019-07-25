function S = pararea(A, B)
% retrun Area of parallelogram(s)

S = sqrt((A(:,2).*B(:,3)-A(:,3).*B(:,2)).^2 + ...
         (A(:,3).*B(:,1)-A(:,1).*B(:,3)).^2 + ...
         (A(:,1).*B(:,2)-A(:,2).*B(:,1)).^2);
return