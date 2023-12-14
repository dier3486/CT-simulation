function r = isavail(x)
% check if x is not empty, not nan, not inf

if isempty(x)
    r = false;
elseif isstruct(x)
    r = ~isempty(fieldnames(x));
else
    r = isfinite(x);
end

end