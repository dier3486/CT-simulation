function [h, map] = objectvisual(object, hfig, color)
% to plot(mesh) an object
% [h, map] = objectvisual(object, hfig, color);
% or h = objectvisual(object);
% e.g.
%   f1=figure; hold on;
%   [h1, map1]=objectvisual(object1,f1,color1);  CD1=h1.CData;
%   [h2, map2]=objectvisual(object2,f1,color2);  CD2=h2.CData;
%   set(h1, 'CData', (CD1+1).*0.9999);  set(h2, 'CData', (CD2+1).*0.9999+2);
%   colormap([map1; map2]);
%   caxis([0,4]);
% to plot two objects in different color 


if nargin<2
    figure();
elseif isempty(hfig)
    figure(gcf);
else
    figure(hfig);
end
if nargin<3 || isempty(color)
    color = hsv2rgb(rand(1, 3).*0.8 + 0.2);
end

% object meshgrid
[X, Y, Z, C] = objectmeshgrid(object);

cscale = linspace(0.6, 1.3, 64);
map = 1-cscale(:)*(1-color);
map(map>1) = 1;
map(map<0) = 0;

colormap(map);
h = mesh(X, Y, Z, C, 'FaceAlpha', 0.3);

end

