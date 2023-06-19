function u = almostlsq(X, Y, b, s)
% an lsq polyfit of binary polynomial, min_u|X'*U*Y-b|^2.
% u = almostlsq(X, Y, b, s);
% INPUT
%       X       X=[x.^0; x.^1; x.^2; ...; x.^(m-1)], x is the x samples
%       Y       Y=[y.^0; y.^1; y.^2; ...; y.^(n-1)], y is the y samples
%       b       b is the target samples
%       s       s is a logical matrix in size (m,n), to select the coefficient to fit
% OUTPUT
%       u       u is the coefficients.
% for U(s)=u, U(~s)=0, the polynomial is f(x,y) = X'*U*Y.

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

s = logical(s);
[m, n] = size(s);
[indexX, indexY] = ndgrid(1:m, 1:n);

% Zp is the outer product of X and Y (s selected)
Zp = X(indexX(s(:)), :).*Y(indexY(s(:)), :);
% inner products on sample space
B = Zp*b(:);    
A = Zp*Zp';

% Au = B;
u = A\B;

end
