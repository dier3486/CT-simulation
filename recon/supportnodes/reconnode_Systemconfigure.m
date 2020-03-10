function [dataflow, prmflow, status] = reconnode_Systemconfigure(dataflow, prmflow, status)
% support node, call system configure from CTsimulation
% [dataflow, prmflow, status] = reconnode_Systemconfigure(dataflow, prmflow, status);

% parameters set in pipe
sysprm = prmflow.pipe.(status.nodename);
if isfield(sysprm, 'cfgfile')
    syscfgfile = sysprm.cfgfile;
else
    % do nothing
    return;
end

% system configure
configure.system = readcfgfile(syscfgfile);
configure = configureclean(configure);
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% load protocol
SYS.protocol = protocolrecon2simu(prmflow.protocol);
SYS = loadprotocol(SYS);

% return
prmflow.SYS = SYS;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function protocol = protocolrecon2simu(protocol)
% we have different defination of these tags:

if isnumeric(protocol.focalspot)
    protocol.focalspot = find(fliplr(dec2bin(protocol.focalspot)=='1'));
end

end