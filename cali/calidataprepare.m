function datafile = calidataprepare(toloop, filepath, fileext, newdatafield)
% prepare data for calibration
% datafile = calidataprepare(toloop, filepath, fileext);
% or, datafile = calidataprepare(toloop, filepath, fileext, newdatafield);
% sorry for no more helps, plz look up the cali scripts for more information.

if nargin<3
    % fileext
    fileext = '.raw';
end
if nargin<4
    % default fieldname
    newdatafield = 'filename';
end

% expand the values in toloop to datafile
datafile = loopcurse(toloop);

% get the file names
datafile = calicorrprepare(datafile, filepath, fileext, newdatafield);

end


function [datafile, idata] = loopcurse(toloop, datafile, ifield, idata)
% expand the values in toloop to datafile

if nargin<2
    datafile = struct();
end
if nargin<3
    ifield = 1;
    idata = 1;
end

toloopfields = fieldnames(toloop);
Nfields = size(toloopfields(:), 1);

currfield = toloop.(toloopfields{ifield});
Nlist = size(currfield(:), 1);
if Nlist<=1 && ~iscell(currfield)
    currfield = {currfield};
    Nlist = 1;
end
% loop the elements in current field
for ii = 1:Nlist
    % copy the values in current field to datafile
    if iscell(currfield)
        datafile(idata).(toloopfields{ifield}) = currfield{ii};
    else
        datafile(idata).(toloopfields{ifield}) = currfield(ii);
    end
    % to recurse
    if ifield<Nfields
        [datafile, idata] = loopcurse(toloop, datafile, ifield+1, idata);
    end
    % next
    if ii<Nlist
        datafile(idata+1) = datafile(idata);
        idata = idata+1;
    end
end

end