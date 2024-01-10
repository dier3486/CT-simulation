function [Zeff, Weff] = effectZandMass(material)
% effective 'atomic' Z and mass of a molecular
% [Zeff, Weff] = effectZandMass(material).
% The Z of a molecular is the sum of its elements, and the mass of it is also the sum of its elements.
% This molecular can be nearly equal to N effective atoms with effective Z=Zeff and effective mass=Weff.

% the molecular is a 'big atom' in Z=sumZ, mass=sumW
sumZ = [material.elemprm(:).Z]*material.elemmol;
sumW = [material.elemprm(:).weight]*material.elemmol;

% meanly devide the 'big atom' to small effective atoms in this way,
Zeff = (([material.elemprm(:).Z].^4*material.elemmol) / sumZ).^(1/3);
Weff = sumW/sumZ * Zeff;

end