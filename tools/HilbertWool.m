function idx = HilbertWool(imagesizeX, imagesizeY)

if nargin == 1
    m0 = [imagesizeX(2); imagesizeX(1)];
else
    m0 = [imagesizeY; imagesizeX];
end
idx0 = 1 : prod(m0);

idx = cummot(idx0, m0);

end

function [idx1, m1] = cummot(idx0, m0)

idx1 = [];
m1 = [];
n = 0;
for ii = 1: size(m0,2)
    m = prod(m0(:, ii));
    [idx_ii, m_ii] = mot(idx0(n+1:n+m), m0(:, ii));
    idx1 = [idx1 idx_ii];
    m1 = [m1 m_ii];
    n = n + m;
end

if any(m1(:)>1)
    [idx1, m1] = cummot(idx1, m1);
end

end

function [idx1, m1] = mot(idx0, m0)

m1 = floor(m0./2);
idx0 = reshape(idx0, m0(:)');

x1 = (idx0(1 : m1(1), 1 : m1(2)))';
x2 = idx0(m1(1)+1 : m0(1), 1 : m1(2));
x3 = idx0(m1(1)+1 : m0(1), m1(2)+1 : m0(2));
x4 = flipud(flipud(idx0(1 : m1(1), m1(2)+1 : m0(2)))');

idx1 = [x1(:); x2(:); x3(:); x4(:)]';

m1 = [m1(2)  m0(1)-m1(1)  m0(1)-m1(1)  m0(2)-m1(2); m1(1)  m1(2)  m0(2)-m1(2)  m1(1)];
m1 = m1(:, all(m1));

end