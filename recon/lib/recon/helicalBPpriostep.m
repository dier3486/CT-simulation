function [dataflow, prmflow, status] = helicalBPpriostep(dataflow, prmflow, status)
% Helical BP priostep

if exist('helicalBPpriostep_advance', 'file')
    % the advaced helicalBPpriostep
    [dataflow, prmflow, status] = helicalBPpriostep_advance(dataflow, prmflow, status);
else
    [dataflow, prmflow, status] = helicalBPpriostep_dump(dataflow, prmflow, status);
end