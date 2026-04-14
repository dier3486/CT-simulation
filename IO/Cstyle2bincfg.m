function bincfg = Cstyle2bincfg(bincfg)

if regexp(class(bincfg), '^lib\.')
    % bincfg is a lib.struct
    bincfg = get(bincfg);
end

bincfg = renameStructField2(bincfg, 'dataclass', 'class');
bincfg = getCchar32(bincfg);

end

function A  = getCchar32(A)

for ifield = fieldnames(A)'
    if strcmp(ifield{1}, 'class')
        C = [A.(ifield{1}); zeros(1, 32)];
        C = typecast(C(:)', 'char');
        C = C(1 : find(C==0, 1, 'first')-1);
        A.(ifield{1}) = C;
    elseif isstruct(A.(ifield{1}))
        A.(ifield{1}) = getCchar32(A.(ifield{1}));
    end
end



end