function x = pathclbs(x)
% clean the back-slashes in path x
% x = pathclbs(x)

x = regexprep(regexprep(x, '\\', '/'), '//', '/');
