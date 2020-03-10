function [dataflow, prmflow, status] = reconnode_simulation(dataflow, prmflow, status)
% support node, call CTsimulation
% [dataflow, prmflow, status] = reconnode_simulation(dataflow, prmflow, status);

% parameters set in pipe
simuprm = prmflow.pipe.(status.nodename);
if isempty(simuprm)
    return;
end

% SYS
SYS = prmflow.SYS;

% simulation pipe line
simupipe = fieldnames(simuprm);
Npipe = size(simupipe(:), 1);
for ii = 1:Npipe
    simunode = simupipe{ii};
    switch lower(simunode)
        case {'projectionscan', 'projection'}
            if isfield(simuprm.(simunode), 'method')
                method = simuprm.(simunode).method;
            else
                method = [];
            end
            Data = projectionscan(SYS, method, 0);
        case 'airprojection'
            if isfield(simuprm.(simunode), 'method')
                method = simuprm.(simunode).method;
            else
                method = [];
            end
            Data = airprojection(SYS, method);
        case 'photon2electron'
            Data = photon2electron(SYS, Data, 0);
        case 'simuresultsoutput'
            simuresultsoutput(SYS, Data);
        otherwise
            error('Unknown simulation function %s', simunode);
    end
end

% return
dataflow.simudata = Data;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
