function [currpool, nextpool] = movepointsaftercopy(currpool, nextpool, readnumber, writenumber, newAvail)

if nargin < 4
    writenumber = readnumber;
end
if nargin < 5
    newAvail = 0;
end

if ~isempty(currpool)
    % move currpool's ReadPoint
    for ii = 1:length(readnumber)
        currpool(ii).ReadPoint = currpool(ii).ReadPoint + readnumber(ii);
        if isfield(currpool(ii), 'ReadViewindex')
            currpool(ii).ReadViewindex = currpool(ii).ReadViewindex + readnumber(ii);
        end
    end
end

if ~isempty(nextpool)
    % move next pool's WritePoint
    for ii = 1:length(writenumber)
        % recode previous WritePoint
        nextpool(ii).prevWritePoint = nextpool(ii).WritePoint;
        % move the WritePoint
        nextpool(ii).WritePoint = nextpool(ii).WritePoint + writenumber(ii);
        WriteClosing = nextpool(ii).WritePoint - nextpool(ii).WriteEnd;
        if isfield(nextpool(ii), 'WriteStuck')
            if writenumber(ii)~=0 && WriteClosing==1
                nextpool(ii).WriteStuck = true;
            end
        end
    end
    % move next pool's AvailPoint
    for ii = 1:length(newAvail)
        nextpool(ii).AvailPoint = nextpool(ii).AvailPoint + newAvail(ii);
    end
end

end