function [currpool, currdata, nextpool, nextdata, jobflag] = ...
    clearoutputarray(currpool, currdata, nextpool, nextdata, viewcommon, echo_onoff)
% pipe-line function to treat the data stucked in waiting array in current pool to be copied to next pool
% [currpool, currdata, nextpool, nextdata, jobflag] = clearoutputarray(currpool, currdata, nextpool, nextdata, viewcommon);
% 
% jobflag = 0:      nothing to copy
% jobflag = 1:      copied all available data
% jobflag = 2:      not all available data be copied (next pool is filled, we should end current node)

if nargin<6
    echo_onoff = false;
end

% % AvailNumber
% if isfield(currpool, 'AvailNumber')
%     AvailPoint = currpool.AvailPoint;
% else
%     AvailPoint = currpool.WritePoint - 1;
% end
% % WriteEnd of outpool
% if isfield(currpool, 'WriteEnd')
%     currWriteEnd = min(currpool.WriteEnd, currpool.poolsize);
% else
%     currWriteEnd = currpool.poolsize;
% end
% % isfilled?
% isallavail = AvailPoint == currWriteEnd;

% try to copy data to next pool
[currpool, nextpool, nextdata, Rinfo] = pooldataoutput2(currpool, currdata, nextpool, nextdata, viewcommon);

% move nextpool's AvailNumber
if isfield(nextpool, 'AvailNumber')
    nextpool.AvailNumber = nextpool.AvailNumber + Rinfo.writenum;
end

if echo_onoff
    fprintf('write to next %d views, ', Rinfo.writenum);
end

% recycle the currpool
[currpool, currdata] = poolrecycle(currpool, currdata);

% % to recycle private buffer while,
% switch currpool.recylestrategy
%     case 0
%         % never recycle
%         % do nothing
%     case 1
%         % always recycle
%         % recycle current pool
%         [currpool, currdata] = poolrecycle2(currpool, currdata, currpool.ReadPoint-1);
%     case 2
%         if currpool.WritePoint >= currpool.warningstage
%             [currpool, currdata] = poolrecycle2(currpool, currdata, currpool.ReadPoint-1);
%         end
%     otherwise
%         error('Illeagal recylestrategy %d!', currpool.recylestrategy);
% end
% 
% % auto unlocked
% WriteStuck = isfield(currpool, 'WriteStuck') && currpool.WriteStuck;
% if WriteStuck && (Rinfo.writenum==Rinfo.datasize_towrite)
%     currpool.WriteStuck = false;
% end

% return jobflag
jobflag = Rinfo.jobflag;

end

