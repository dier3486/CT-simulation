% test script 

addpath(genpath('../'));

A = randn(5,3);
B = randn(5,3);

object1.type = 'sphere';
% object1.O = randn(1,3)./10;
% object1.Vector = rand(3);
object1.O = [0, 0, 0];
object1.Vector = eye(3);
object1.invV = inv(object1.Vector);
object1.ref_mu = 0.01;
object1.material = 'water';

[D, L] = intersection(A, B, object1);