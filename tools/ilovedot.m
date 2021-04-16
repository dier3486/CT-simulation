function str_dot = ilovedot(str)
% add '.' before '*','/' and '^' in str

str_dot = regexprep(str, '(?=[^\.])[\*\/\^]', '\.$0');

end