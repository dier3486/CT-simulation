function strlinked = linkstringcell(strcell, linkkey)
% to link a cell of strings
% strlinked = linkstringcell(strcell, linkkey);
% or 
% strlinked = linkstringcell(strcell);
% e.g. linkstringcell({'a', 'bc'}, '_') = '_a_bc';

if nargin<2
    linkkey = '';
end

if isempty(strcell)
    strlinked = '';
    return
end

strcell(2, :) = strcell;
strcell(1, :) = {linkkey};
strlinked = [strcell{:}];

end

