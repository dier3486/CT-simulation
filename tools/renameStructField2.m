function STR = renameStructField2(STR, OLDFIELDNAME, NEWFIELDNAME)
% recursed renameStructField

% call renameStructField
STR = renameStructField(STR, OLDFIELDNAME, NEWFIELDNAME);

% recurse
for ifield = fieldnames(STR)'
    if isstruct(STR.(ifield{1}))
        STR.(ifield{1}) = renameStructField2(STR.(ifield{1}), OLDFIELDNAME, NEWFIELDNAME);
    end
end

end