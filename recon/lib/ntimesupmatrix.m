function A = ntimesupmatrix(m, n, c)
% n-times omiga4 upsampling matrix, with repeat boundary condition
%   A = ntimesupmatrix(m, n, [c1 c2]);
% m is the data length, n is the upsampling times, [c1 c2] the gamma-parameter.
% Then the upsampling result reads,
%   y = x*A;
% x is a matrix in (Nx, m) and y is in (Nx, (m-1)*n+1).
% Note that the A is a full matrix therefore the m shall be a small number. 
% Typically the m is the slice number and Nx is the pixel number to do an upsampling on slice direction.

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

A = zeros(m, (m-1)*n+1, 'like', c);

n = cast(n, 'like', c);
alpha = (0:n)'./n;
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
gamma = c(1)./sqrt(1-alpha.*(1-alpha).*c(2));
t = [(2-alpha-beta+(2-alpha.*3).*gamma)./4 (beta+(1-alpha).*(1-gamma))./4];

for ii = 1:n
    % 1
    for jj = 0 : (m - 1 - (ii~=1))
        p = jj + (jj==0);
        q = (jj*n) + ii;
        A(p, q) = A(p, q) + t(ii, 2);
    end
    % 2
    for jj = 1 : (m - (ii~=1))
        p = jj;
        q = (jj-1)*n + ii;
        A(p, q) = A(p, q) + t(ii, 1);
    end
    % 3
    for jj = 2 : (m + (ii==1))
        p = jj - (jj>m);
        q = (jj-2)*n + ii;
        A(p, q) = A(p, q) + t(n+2-ii, 1);
    end
    % 4
    for jj = 3 : m+1
        p = jj - (jj>m);
        q = (jj-3)*n + ii;
        A(p, q) = A(p, q) + t(n+2-ii, 2);
    end

end

end