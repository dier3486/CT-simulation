function jstruct = jsonread(fname)
% read json file

f1 = fopen(fname, 'r');
charf1 = fread(f1, inf, 'char=>char');
fclose(f1);
jstruct = jsondecode(charf1(:)');
