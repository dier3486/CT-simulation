function cs = pinsplines(rawdata, fanangles, Npixel, Nslice, Nview)
% % subroutine used in geometric cali 

rawdata = reshape(rawdata, Npixel, []);
[pmax, idxmax] = max(rawdata, [], 1);

m = 20;
index = (-m:m)' + idxmax;
sout = (index>Npixel) | (index<1);
index(index>Npixel) = Npixel;
index(index<1) = 1;
index = index + (0:Nview*Nslice-1).*Npixel;
A1 = rawdata(index);
A1(sout) = 0;
C1 = (sum(rawdata, 1) - sum(A1, 1))./(Npixel-m*2-1-sum(sout, 1));

A1 = A1 - C1;
A1(sout) = 0;
pmax = pmax - C1;

cut = 0.05;
A1(A1 < repmat(pmax.*cut, m*2+1, 1)) = 0;

cs = spline([0 1], [0 0]);
cs.vindex = [];
cs.pindex = [];
cs.p = nan;
cs.dp = nan;
cs(Nslice, Nview) = cs(1);
for ii = 1:Nview*Nslice
    islice = mod(ii-1, Nslice) + 1;
    iview = ceil(ii/Nslice);
    find_ii = find(A1(:, ii)>0, 1, 'first')-1 : find(A1(:, ii)>0, 1, 'last')+1;
    fanindex_ii = mod(index(find_ii, ii)'-1, Npixel)+1;
    fanindex_ii = unique(fanindex_ii);
    fan_ii = fanangles(fanindex_ii, islice);
    if length(fanindex_ii) == length(find_ii)
        tmp = spline(fan_ii, [0; A1(find_ii, ii); 0]);
        tmp.pindex = fanindex_ii;
        tmp.vindex = ones(size(fanindex_ii)).*iview;
        [tmp.dp, tmp.p] = ppcenterderivative(tmp);
        cs(ii) = tmp;
    else
        % pin out of FOV
        cs(ii).pindex = fanindex_ii;
        cs(ii).vindex = ones(size(fanindex_ii)).*iview;
        cs(ii).p = nan;
        cs(ii).dp = zeros(size(fanindex_ii));
    end
end

end