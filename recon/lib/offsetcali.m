function offset = offsetcali(offset, rawdata, w)

if isempty(offset)
    offset = struct();
    offset.Nview = 0;
    offset.variance = 0;
    offset.main = 0;
end

if nargin < 3
    % var weight
    w = true;
end

Nviewpre = offset.Nview;
Nviewcurr = size(rawdata, 2);

maincurr = mean(rawdata, 2);
varcurr = var(rawdata, w, 2);

nw = ~w;
offset.variance = offset.variance .* ((Nviewpre-nw)/(Nviewpre + Nviewcurr-nw)) +...
                  varcurr .* ((Nviewcurr-nw)/(Nviewpre + Nviewcurr-nw)) + ...
                  (offset.main-maincurr).^2 .* (Nviewpre*Nviewcurr/(Nviewpre + Nviewcurr-nw)/(Nviewpre + Nviewcurr));
offset.main = offset.main + (maincurr - offset.main) .* (Nviewcurr/(Nviewcurr + Nviewpre));
offset.Nview = Nviewpre + Nviewcurr;
end