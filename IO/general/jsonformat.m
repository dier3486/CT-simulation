function [jschar, idx] = jsonformat(jschar, idx, bracelv)
% add tab, blank and enters to jschar for easy to read
% only support struct frame (start with '{') but not cell frame (start with
% '['.

if nargin < 2
    idx = 1;
end
if nargin < 3
    bracelv = 0;
end
ent = newline();
tab = sprintf('\t');
inbracket = 0;
inquot = false;
slashed = false;
while idx<=length(jschar)
    if jschar(idx) == '\'
        slashed = ~slashed;
    else
        slashed = false;
    end
    switch jschar(idx)
        case '"'
            if ~slashed
                inquot = ~inquot;
            end
        case '['
            if ~inquot
                inbracket = inbracket+1;
            end
        case ']'
            if ~inquot
                inbracket = inbracket-1;
            end
        case '{'
            if ~inquot
                tabs = repmat(tab, 1, bracelv);
                inserts = ['{', ent, tabs];
                jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
                idx = idx + bracelv + 1;
                % recurse
                [jschar, idx] = jsonformat(jschar, idx, bracelv + 1);
            end
        case '}'
            if ~inquot
                tabs = repmat(tab, 1, bracelv - 1);
                inserts = [ent, tabs, '}'];
                jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
                idx = idx + bracelv;
                return
            end
        case ':'
            if ~inquot
                inserts = ': ';
                jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
                idx = idx + 1;
            end
        case ','
            if ~inquot
                if inbracket>0
                    inserts = ', ';
                    jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
                    idx = idx + 1;
                else
                    tabs = repmat(tab, 1, bracelv - 1);
                    inserts = [',', ent, tabs];
                    jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
                    idx = idx + bracelv;
                end
            end
        otherwise
            1;
    end
    idx = idx + 1;
end