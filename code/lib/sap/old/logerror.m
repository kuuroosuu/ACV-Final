function logerror = logerror(Data, index, logfile, reason)
% LOGERROR Log error midi/wave files

if isempty(index), return; end

fid = fopen(logfile, 'a');
fprintf(fid, '%s\n', ['====== (', int2str(length(index)), ') ', reason, ':']);
if isfield(Data, 'midiName'),	% songData
	for i=1:length(index),
		fprintf(fid, '%d. %s\n', index(i), Data(index(i)).midiName);
	end
else	% waveData
	for i=1:length(index),
		fprintf(fid, '%d. %s ===> %s (%s)\n', index(i), Data(index(i)).name, Data(index(i)).errorMsg, reason);
	end
end
fprintf(fid, '\n');
fclose(fid);