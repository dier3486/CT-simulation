function jschar = jsonwrite(jstruct, fname)
% write json file
% jsonwrite(jstruct, fname)

jschar = jsonformat(jsonencode(jstruct));
f1 = fopen(fname, 'w');
if f1>=0
    fwrite(f1, jschar);
end
if f1>=3
    fclose(f1);
end

return
