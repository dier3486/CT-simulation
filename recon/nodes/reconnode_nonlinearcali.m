function [dataflow, prmflow, status] = reconnode_nonlinearcali(dataflow, prmflow, status)
% nonlinear calibration
% [dataflow, prmflow, status] = reconnode_nonlinearcali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;

% parameters to use
if ~isempty(status)
    caliprm = prmflow.pipe.(status.nodename);
else
    % for debug
    caliprm = struct();     
end

% polynomial order
if isfield(caliprm, 'polyorder')
    n_ploy = caliprm.polyorder;
else
    % order of nonlinear poly, default is 2
    n_ploy = 2;
    % NOTE: large order could be unstable
end
% HCscale
if isfield(caliprm, 'HCscale')
    HCscale = caliprm.HCscale;
else
    HCscale = 1000;
end
% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end
% echo .
echo_onoff = true;

% I know the input data are rawdata_bk1, rawdata_bk2 ...
% find out them
datafields = findfields(dataflow, '\<rawdata_bk');
Nbk = length(datafields);
% I know Nbk must be even
% and the rawhead
headfields = findfields(dataflow, '\<rawhead_bk');

% reshape
for ii = 1:Nbk
    dataflow.(datafields{ii}) = reshape(dataflow.(datafields{ii}), Nps, Nview);
end

% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

% ini the buffer to use
x = zeros(Nbk/2, Nview);
y = zeros(Nbk/2, Nview);
% ini the coefficients of the polynomal
poly_nonl = nan(Npixel*Nslice, n_ploy);
% fitting option
options = optimoptions('lsqnonlin','Display','off');
% initial value of the fitting
t0 = zeros(1, n_ploy);
t0(n_ploy) = 1.0;

% loop the pixels to fit the polynomal for each
for islice = 1:Nslice
    if echo_onoff, fprintf('.'); end
    % index range
    Srange = zeros(Npixel, Nview, Nbk/2);
    for ibk = 1:Nbk/2
        for iview = 1:Nview
            Srange(index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk), iview, ibk) = 1;
        end
    end
    
    % For each slice we will start from the middle pixel that will be more stable though a little troublesome.
    midu = floor(Npixel/2);
    % ini the initial value for ipixel
    t0_ipx = t0;
    % loop from midu+1 to end, then from midu to 1
    loopindex = [(midu+1:Npixel) (midu:-1:1)];
%     % debug 
%     loopindex = [(midu+1:midu+5) (midu:-1:midu-5)];
    for ii = loopindex
        % reset t0 when restart from midu
        if ii == midu
            t0_ipx = t0;
        end
        % pixel index
        ipixel = ii+(islice-1)*Npixel;
        % set the data to fit
        for ibk = 1:Nbk/2
            % x is the original data in fitting
            ix = ibk*2-1;   % 1 3 5 7
            x(ibk, :) = double(dataflow.(datafields{ix})(ipixel, :)) ./ HCscale;
            % y is the target data
            iy = ibk*2;     % 2 4 6 8
            y(ibk, :) = double(dataflow.(datafields{iy})(ipixel, :)) ./ HCscale;
            % cut y by index range
            y(ibk, :) = y(ibk, :).*Srange(ii, :, ibk);
        end
        % Is there enough available data for fitting?
        s = y>0;
        s_sum = sum(s, 2);
        if all(s_sum<Nview/4)
            % skip this pixel 
            continue;
        end
        % to find out a function f: min(|f(x)-y|)
        % the algorithm is reverse generator base on lsqnonlin
        p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)), t0_ipx, [], [], options);
        if p_ipx(1)>0.5
            % unstable? redo the fitting start with t0
            p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)), t0, [], [], options);
        end
        % set the initial value for next pixel
        t0_ipx = p_ipx;
        % copy to result buffer
        poly_nonl(ipixel, :) = p_ipx;
    end
    
end
% fillup nan
poly_nonl = reshape(fillmissing(reshape(poly_nonl, Npixel, []), 'nearest'), [], n_ploy);

% paramters for corr
nonlinearcorr = caliprmforcorr(prmflow, corrversion);
% copy results to corr
nonlinearcorr.Nslice = Nslice;
nonlinearcorr.order = n_ploy;
nonlinearcorr.mainsize = Nps*n_ploy;
nonlinearcorr.main = poly_nonl;

% to return
dataflow.nonlinearcorr = nonlinearcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
