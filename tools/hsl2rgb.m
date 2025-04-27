function C = hsl2rgb(C)
% HSL to RGB

C = hsv2rgb(hsl2hsv(C));

end
