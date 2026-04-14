function [u, b, delta] = BregmanTV(img0, mu, Cl, u0, b, varargin)
% TV

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

TVopts = TVoptions(varargin{:});

Ndim = str2double(regexp(TVopts.DIM, '^\d', 'match'));
balance = TVopts.Balance;
% if length(balance) < Ndim
%     balance(end+1 : Ndim) = 1;
% else
%     balance = balance(1 : Ndim);
% end
% balance = balance.*(Ndim/sum(balance));
% if length(balance) < 3
%     balance(end+1 : 3) = 0;
% end
% boundary condition
BC = TVboundarycond(TVopts);

fclass = classGPU(img0);

% f
f1 = img0;
Srange = (f1>=TVopts.Crange(1)) & (f1<=TVopts.Crange(2));
Nsrang = sum(Srange(:));
f1(f1<TVopts.Crange(1)) = TVopts.Crange(1);
f1(f1>TVopts.Crange(2)) = TVopts.Crange(2);

imgsize = size(img0);
if length(imgsize) < 3
    imgsize(end+1 : 3) = 1;
end

% initial u0
if isempty(u0)
    u0 = f1;
end

% initial d b
if isempty(b)
    b = zeros([imgsize 3], 'like', img0);
    d = zeros([imgsize 3], 'like', img0);
else
    [d, b] = TVfundbyub(u0, b, mu.*Cl, Ndim, BC);
end

% iteration
delta = zeros(1, TVopts.MaxIter, fclass);
for ii = 1 : TVopts.MaxIter
    % u
    if ii > 1
        u0 = u;
    end
    % funG
    if exist('TVfunGadv', 'file') == 2
        u = TVfunGadv(f1, u0, b, d, mu, Cl, balance, Ndim, TVopts, BC);
    else
        % the advanced TV version may be not authorized
        u = TVfunG(f1, u0, b, d, mu, Cl, balance, Ndim, TVopts, BC);
    end
    % delta
    delta(ii) = gather(sqrt(sum((u(:)-u0(:)).^2)./Nsrang));
    if delta(ii) < TVopts.tol || ii == TVopts.MaxIter
        break;
    end
    % d b
    [d, b] = TVfundbyub(u, b, mu.*Cl, Ndim, BC);
end

if TVopts.MaxIter>1 && delta(end)>delta(end-1)
    warning('TV failed!');
end

u(~Srange) = img0(~Srange);

end
