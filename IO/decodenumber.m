function r = decodenumber(c, Sorig, Sprev)
% a subfunction to explain the numers in cfg_ii.size and cfg_ii.number used in sparsepack.m and packstruct.m

if nargin < 3
    Sprev = Sorig;
end

if isnumeric(c)
    r = c;
elseif isempty(c)
    r = [];
elseif ischar(c)
    %         c(c=='$') = 'S';
    c = regexprep(c, '\$\$' , 'Sorig');
    c = regexprep(c, '\$' , 'Sprev');
    try
        r = eval(c);
    catch
        r = c;
    end
else
    r = nan;
end
% to double
if isnumeric(r)
    r = double(r);
end

end