function [currpool, nextpool, nextdata, Rinfo] = pooldataoutput2(currpool, currdata, nextpool, nextdata, ...
    viewcommon)
% to copy the data from currpool to nextpool, inner version
% a simplified version in using pooldatacopy
% [currpool, currdata, nextpool, nextdata, output] = pooldataoutput2(currpool, currdata, nextpool, nextdata);

% the view number must be integer times to viewcommon
if nargin < 5
    viewcommon = 1;
end

% current pool circulatation mode
if isfield(currpool, 'circulatemode')
    currcirculate = currpool.circulatemode;
else
    currcirculate = false;
end

% next pool circulatation mode
if isfield(nextpool, 'circulatemode')
    nextcirculate = nextpool.circulatemode;
else
    nextcirculate = false;
end

% check error in currpool.AvailNumber
if currcirculate
    if currpool.AvailNumber > currpool.poolsize
        currpool.AvailNumber = currpool.poolsize;
        warning('Due to the AvailNumber is out of the poolsize it is forced to set to the most limitation!');
    end
else
    if currpool.AvailNumber > currpool.poolsize - currpool.ReadPoint + 1
        currpool.AvailNumber = currpool.poolsize - currpool.ReadPoint + 1;
        warning('Due to the AvailNumber is out of the poolsize it is forced to set to the most limitation!');
    end
end

% datasize to write
RP2RE = poolp2p(currpool.ReadPoint, currpool.ReadEnd, currpool.poolsize, currcirculate);
datasize_towrite = min(currpool.AvailNumber, RP2RE);
datasize_towrite = max(datasize_towrite, 0);
datasize_towrite = floor(datasize_towrite/viewcommon) * viewcommon;

% % shot limit
% if shotlimit && isfield(currdata, 'rawhead')
%     % view index to read
%     readingindex = currpool.ReadPoint : currpool.ReadPoint + datasize_towrite-1;
%     if currcirculate
%         readingindex = mod(readingindex - 1, currpool.poolsize) + 1;
%     end
%     % shot index to read
%     Shot_Number = currdata.rawhead.Shot_Number(readingindex);
%     end_shot = find(Shot_Number ~= Shot_Number(1), 1, 'first');
%     if ~isempty(s_shot)
%         datasize_towrite = end_shot - 1;
%     end
% end
% We will use a circulating pool and writelock flag to stuck next shot data
% untill the current shot is cleared; for multi shots helical or half axial
% or other non-axial multi shots cases, the stuck function will be employed
% by the pool.ReadEnd but not to check the rawhead.
% So we abandoned the 'shot limit' function here. No any rawdata depending 
% codes shall be employed in this function. 

% datasize_towrite will be returned in Rinfo
Rinfo.datasize_towrite = datasize_towrite;

if ~isempty(nextpool)
    % the original WritePoint of next pool will be Rinfo
    Rinfo.WritePoint = nextpool.WritePoint;
    % is writing forbidened
    WriteStuck = isfield(nextpool, 'WriteStuck') && nextpool.WriteStuck;
else
    Rinfo.WritePoint = 0;
    WriteStuck = false;
end

% is reading forbidened
ReadStuck =  isfield(currpool, 'ReadStuck') && currpool.ReadStuck;

% stucked
if WriteStuck || ReadStuck
    Rinfo.writenum = 0;
    return;
end

if currcirculate && ...
        (datasize_towrite > currpool.poolsize - currpool.ReadPoint + 1)
    % circulately out of size
    % record AvailNumber
    recAvailNumber = currpool.AvailNumber;
    % reset AvailNumber to 
    currpool.AvailNumber = currpool.poolsize - currpool.ReadPoint + 1;
    % recurse
    [currpool, nextpool, nextdata, Rinfo_rec] = ...
        pooldataoutput2(currpool, currdata, nextpool, nextdata, viewcommon);
    % recover AvailNumber
    currpool.AvailNumber = recAvailNumber - Rinfo_rec.writenum;
    % new datasize_towrite
    datasize_towrite = datasize_towrite - Rinfo_rec.writenum;
    % prepare to return the writenum
    Rinfo.writenum = Rinfo_rec.writenum;
else
    Rinfo.writenum = 0;
end

if ~isempty(nextpool)
    % WritePoint to WriteEnd
    WP2WE = poolp2p(nextpool.WritePoint, nextpool.WriteEnd, nextpool.poolsize, nextcirculate);
    % writenum
    writenum = min([datasize_towrite, nextpool.poolsize - nextpool.WritePoint + 1, WP2WE+1]);
    writenum = max(writenum, 0);
    writenum = floor(writenum/viewcommon) * viewcommon;

    % accept fields
    if isfield(nextpool, 'datafields')
        datafields = nextpool.datafields;
    else
        % default is all
        datafields = {};
    end
    1;
    % copy curr pool to next pool
    [nextdata, writenum] = ...
        pooldatacopy(currdata, nextdata, currpool.ReadPoint, nextpool.WritePoint, writenum, datafields, true);
    
    % move next pool's write point
    nextpool.WritePoint = nextpool.WritePoint + writenum;
    if nextcirculate
        nextpool.WritePoint = mod(nextpool.WritePoint-1, nextpool.poolsize) + 1;
    end

else
    % Even when the nextnode is not existing the node will do its work as usual.
    writenum = datasize_towrite;
    Rinfo.WritePoint = 0;
end

% move currpool's read point
currpool.ReadPoint = currpool.ReadPoint + writenum;
if currcirculate
    currpool.ReadPoint = mod(currpool.ReadPoint-1, currpool.poolsize) + 1;
end
currpool.AvailNumber = currpool.AvailNumber - writenum;
if isfield(currpool, 'ReadViewindex')
    currpool.ReadViewindex = currpool.ReadViewindex + writenum;
end

% ReadEnd to ReadPoint of currpool
ReadingClosing = poolp2p(currpool.ReadEnd, currpool.ReadPoint, currpool.poolsize, currcirculate);
% WriteEnd to WriteEnd of nextpool
WriteClosing = poolp2p(nextpool.WriteEnd, nextpool.WritePoint, nextpool.poolsize, nextcirculate);

% read/write self lock
if currcirculate && isfield(currpool, 'ReadStuck')
    currpool.ReadStuck = (writenum>0) & (ReadingClosing==1);
end
if nextcirculate && isfield(nextpool, 'WriteStuck')
    nextpool.WriteStuck = (writenum>0) & (WriteClosing==1);
end

% when next pool is circulatation try to recurse
if nextcirculate && (writenum<datasize_towrite) && (WriteClosing~=1)
%     nextpool.WritePoint = 1;
    % recurse
    [currpool, nextpool, nextdata, Rinfo_rec] = pooldataoutput2(currpool, currdata, nextpool, nextdata, ...
        viewcommon);
    % do not to recurse twice!
    % new writenum
    writenum = writenum + Rinfo_rec.writenum;
end

% to return writenum
Rinfo.writenum = Rinfo.writenum + writenum;

% flag of jobdone
if (Rinfo.datasize_towrite > 0) && (Rinfo.datasize_towrite > Rinfo.writenum)
    % nextpool could be stucked
    Rinfo.jobflag = 2;
elseif Rinfo.datasize_towrite > 0
    % wrote something
    Rinfo.jobflag = 1;
else
    % nothing wrote
    Rinfo.jobflag = 0;
end

end