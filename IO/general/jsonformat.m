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
ent = sprintf('\n');
tab = sprintf('\t');
inbracket = 0;
while idx<=length(jschar)
    switch jschar(idx)
        case '['
            inbracket = inbracket+1;
        case ']'
            inbracket = inbracket-1;
        case '{'
            tabs = repmat(tab, 1, bracelv);
            inserts = ['{', ent, tabs];
            jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
            idx = idx + bracelv + 1;
            % recurse
            [jschar, idx] = jsonformat(jschar, idx, bracelv + 1);
%             bracelv = bracelv + 1;
        case '}'
            tabs = repmat(tab, 1, bracelv - 1);
            inserts = [ent, tabs, '}'];
            jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
            idx = idx + bracelv;
            return
%             bracelv = bracelv - 1;
        case ':'
            inserts = ': ';
            jschar = [jschar(1:idx-1), inserts, jschar(idx+1:end)];
            idx = idx + 1;
        case ','
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
    idx = idx + 1;
end