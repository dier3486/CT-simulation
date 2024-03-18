function r = splitcomma(string, mycomma)
% regexp split ',' (or mycomma) of a string
% r = splitcomma(string, mycomma);
% or, r = splitcomma(string);
% It is regexp(string, '(\s+,\s+)|(\s+,)|(,\s+)|(,)|(^\s+)|(\s+$)', 'split');

if nargin < 2
    expr = '(\s+,\s+)|(\s+,)|(,\s+)|(,)|(^\s+)|(\s+$)';
else
    expr = sprintf('(\\s+%s\\s+)|(\\s+%s)|(%s\\s+)|(%s)|(^\\s+)|(\\s+$)',mycomma,mycomma,mycomma,mycomma);
end

r = regexp(string, expr, 'split');
s = cellfun(@(x)~isempty(x), r);
r = r(s);

end