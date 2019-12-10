function [prmflow, status] = reconinitial(prmflow, status)
% recon initial
% [prmflow, status] = reconinitial(prmflow, status)

prmflow.system = status.protocol.system;
status.series_index = 1;

return