function [indexIn, indexOut] = indexinoutmap(Nin, Nout, Rely, D)
% help to set the maping index from the input to output 

indexIn = max(1 + Rely(2) + D, 1) : min(Nout + Rely(2) + D, Nin + Rely(1) + Rely(2));
indexOut = max(1 - Rely(2) - D, 1) : min(Nout, Nin + Rely(1) - D);


end