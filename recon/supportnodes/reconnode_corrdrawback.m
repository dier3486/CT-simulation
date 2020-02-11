function [dataflow, prmflow, status] = reconnode_corrdrawback(dataflow, prmflow, status)
% recon node, draw back correction(s)
% [dataflow, prmflow, status] = reconnode_corrdrawback(dataflow, prmflow, status);
% used in some calibration, yep it is totally crazy

% parameters set in pipe
drawbackprm = prmflow.pipe.(status.nodename);

dbfields = fieldnames(drawbackprm);
for ii = 1:length(dbfields)
    drawbacknode = dbfields{ii};
    nodename_slip = regexp(drawbacknode, '_', 'split');
    switch lower(nodename_slip{1})
        case 'housefield'
            % drawback Housefield correction
            if isfield(drawbackprm.(drawbacknode), 'HCscale')
                HCscale = drawbackprm.(drawbacknode).HCscale;
            elseif isfield(prmflow.pipe, drawbacknode) && isfield(prmflow.pipe.(drawbacknode), 'HCscale')
                HCscale = prmflow.pipe.(drawbacknode).HCscale;
            else
                HCscale = 1000;
            end
            % drawback the scale
            dataflow.rawdata = dataflow.rawdata./HCscale;
            % NOTE: the housefield node in prflow.pipe could be flexibly named, e.g. prflow.pipe.Housefield_kickass, in this
            % case if we want to cite whose parameters and/or corrtable, which is 'prmflow.corrtable.Housefield_kickass', we
            % shall set the node in corrdrawback with same name, like prflow.pipe.corrdrawback.Housefield_kickass
            % And, the parameters in pipe can be replaced, e.g. to set pipe.corrdrawback.Housefield_kickass.HCscale=1024 to
            % set a new HCscale value, be the corrtable can not be replaced, at least in this function can not. If you really 
            % want to commit that, there are some very helpful virus nodes in this folder..
        case {'beamharden', 'nonlinear'}
            % TBC
            1;
        case 'crosstalk'
            % TBC
            1;
        case {'air', 'aircorr'}
            % TBC
            1;
        case 'log2'
            % draw back log2
            % NOTE: it is not exp(dataflow.rawdata), it is really drawback
            dataflow = drawbacklog2(drawbackprm.(drawbacknode), dataflow, prmflow);
            % 
        case 'filter'
            % TBC
            1;
        case 'axialrebin'
            % plz call reconnode_inverserebin.m instead
            1;
        otherwise
            % not support
            warning('%s is not supported to drawback.', nodename_slip{1});
    end
end


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function dataflow = drawbacklog2(drawbackprm, dataflow, prmflow)
% drawback log2

% Z0
if isfiled(drawbackprm, 'DBBzero')
    Z0 = drawbackprm.DBBzero;
elseif isfield(prmflow, 'system') && isfield(prmflow.system, 'DBBzero')
    Z0 = prmflow.system.DBBzero;
else
    Z0 = 16384;
end

% reshape
Nview = prmflow.recon.Nview;
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);

% inverse log2
dataflow.rawdata = 2.^(-dataflow.rawdata).*single(dataflow.rawhead.Integration_Time);

% offset
if isfield(dataflow, 'offset')
    dataflow.rawdata = dataflow.rawdata + mean(dataflow.offset.rawdata, 2);
else
    dataflow.rawdata = dataflow.rawdata + Z0;
end

end