function timestamp = posixtimestamp(t)

if nargin<1
    t = 'now';
end

timestamp = int64(posixtime(datetime(t)));

end