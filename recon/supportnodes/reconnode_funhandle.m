function [dataflow, prmflow, status] = reconnode_funhandle(dataflow, prmflow, status)
% virus node, run a function handle on dataflow.rawdata (not always)
% [dataflow, prmflow, status] = reconnode_funhandle(dataflow, prmflow, status);
% WARN: which can do anything

% parameters set in pipe
handleprm = prmflow.pipe.(status.nodename);

if ~isfield(handleprm, 'fun')
    % do nothing
    return
else
    myfun = handleprm.fun;
    if ~isa(myfun, 'function_handle')
        myfun = str2func(myfun);
    end
end
if isfield(handleprm, 'argin')
    if iscell(handleprm.argin)
        argin = handleprm.argin;
    else
        argin = {handleprm.argin};
    end
else
    argin = {};
end
if isfield(handleprm, 'argout')
    if iscell(handleprm.argout)
        argout = handleprm.argout;
    else
        argout = {handleprm.argout};
    end
else
    argout = {};
end

% clean inputs
for ii = 1:length(argin)
    try
        argin{ii} = eval(argin{ii});
    catch
    end
end

% call myfun
funinputs = [dataflow.rawdata argin];
if isempty(argout)
    dataflow.rawdata = feval(myfun, funinputs{:});
else
    Nout = length(argout);
    funoutputs = cell(1, Nout);
    funoutputs{:} = feval(myfun, funinputs{:});
    % assign the outputs to variables
    for ii = 1:Nout
        expression = [argout{ii} ' = funoutputs{' num2str(ii) '}'];
        try
            eval(expression);
        catch
        end
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end