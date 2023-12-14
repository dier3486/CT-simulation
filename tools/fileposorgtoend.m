function datasize = fileposorgtoend(fid)
% datasize left in the file

pforigin = ftell(fid);
fseek(fid, 0, 'eof');
pfend = ftell(fid);
fseek(fid, pforigin, 'bof');

datasize = pfend - pforigin;

end