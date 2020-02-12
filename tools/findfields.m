function fd = findfields(A, tofind)
% to find the fields in A match with the regular expression tofind
% fd = findfields(A, tofind);
% e.g. when tofind='\<rawdata_bk', it will return all the fields' name of A.rawdata_bk*.
% the tofind is a regular expression, plz help regexp for more information

% list the fields name in A
Afields = fieldnames(A);
% if they match with tofind?
findbk = regexp(Afields, tofind);
% find out each of them
findex = false(1, length(Afields));
for ii = 1:length(Afields)
    if ~isempty(findbk{ii})
        findex(ii) = true;
    end
end
% so we find out the fields (names)
fd = Afields(findex);

return