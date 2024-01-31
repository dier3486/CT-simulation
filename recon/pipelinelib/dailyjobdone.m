function jobdone = dailyjobdone(datasize_read, datasize_wrote, datasize_done)
% daily job done
% jobdone = dailyjobdone(datasize_read, datasize_wrote, datasize_done)

if nargin < 3  || isempty(datasize_done)
    datasize_done = datasize_wrote;
end

if datasize_done == 0 && datasize_read == 0
    % no reading no writting, did nothing
    % pass, sleep
    jobdone = 3;
elseif datasize_done == 0
    % read something but make nothing done
    % technically pass, to recycle (the input pool) and sleep
    jobdone = 4;
elseif datasize_wrote < datasize_done
    % read something, did somthing but partly written to next pool
    % technically done, to recycle (the input pool) and wake up next node, keep waking
    jobdone = 2;
else
    % read something and write all available data 
    % normally done
    jobdone = 1;
end

end