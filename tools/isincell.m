function r = isincell(C, x)
% is x in cell C

r = any(strcmp(C(:), x));

end