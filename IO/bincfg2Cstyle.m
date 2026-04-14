function bincfg = bincfg2Cstyle(bincfg)

bincfg = renameStructField2(bincfg, 'class', 'dataclass');
bincfg = everything2single(bincfg, 'double', 'int32');
bincfg = char2Cchar32(bincfg);

end

function A  = char2Cchar32(A)

for ifield = fieldnames(A)'
    if strcmp(ifield{1}, 'dataclass')
        % C = [A.(ifield{1}); zeros(1, 32)];
        % C = typecast(C(:)', 'char');
        % C = C(1 : find(C==0, 1, 'first')-1);
        A.(ifield{1}) = int32(0);
    elseif isstruct(A.(ifield{1}))
        A.(ifield{1}) = char2Cchar32(A.(ifield{1}));
    end
end

% if isstruct(A)
%     for ifield = fieldnames(A)'
%         A.(ifield{1}) = char2Cchar32(A.(ifield{1}));
%     end
% elseif ischar(A)
%     A = typecast(A, 'uint8');
%     % I know the char of matlab is in UCS-2.
%     % cast to ASCII
%     A = A(1:2:end);
%     if length(A) > 32
%         A = A(1:32);
%     elseif length(A) < 32
%         A(32) = 0;
%     end
% end

end