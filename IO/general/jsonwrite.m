function jschar = jsonwrite(jstruct, fname)
% write json file
% jschar = jsonwrite(jstruct, fname);

jstruct = jsonclear(jstruct);
jschar = jsonformat(jsonencode(jstruct));

if nargin > 1
    f1 = fopen(fname, 'w');
    if f1>=0
        fwrite(f1, jschar);
    end
    if f1>=3
        fclose(f1);
    end
end

end

function S = jsonclear(S)

for fields = fieldnames(S)'
    if isgpuarray(S.(fields{1}))
        S.(fields{1}) = gather(S.(fields{1}));
    end
    if isstruct(S.(fields{1}))
        for ii = 1:length(S.(fields{1})(:))
            S.(fields{1})(ii) = jsonclear(S.(fields{1})(ii));
        end
    elseif iscell(S.(fields{1}))
        for ii = 1:length(S.(fields{1})(:))
            S.(fields{1}){ii} = jsonclear(S.(fields{1}){ii});
        end
    else
        switch class(S.(fields{1}))
            case 'function_handle'
                S.(fields{1}) = char(S.(fields{1}));
            otherwise     
        end
    end
end

end