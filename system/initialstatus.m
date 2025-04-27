function status = initialstatus(reconcfg)

status = struct();
status.reconcfg = reconcfg;
status.taskUID = dicomuid();

if isfield(status.reconcfg, 'debug')
    status.debug = status.reconcfg.debug;
else
    status.debug = struct();
end

% default modes
% echo
if ~isfield(status.debug, 'echo_onoff')
    status.debug.echo_onoff = true;
end
% trace of the pools
if ~isfield(status.debug, 'pooltrace_onoff')
    status.debug.pooltrace_onoff = true;
end


end