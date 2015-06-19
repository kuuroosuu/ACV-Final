function success=waveFileRecord(duration, fs, nbits, waveFile, text4record, forcedRecording)
% waveFileRecord: Record from mic and save the data to a wave file
%
%	Usage:
%		success=waveFileRecord(duration, fs, nbits, waveFile, text4record, forcedRecording)
%			forcedRecording: Forced recording without skipping
%
%	Example:
%		duration=3;
%		fs=16000;
%		nbits=16;
%		waveFile=[tempname, 'test.wav'];
%		text4record='This is an example text';
%		waveFileRecord(duration, fs, nbits, waveFile, text4record);

%	Roger Jang, 20040424, 20130325, 20150320

if nargin<1, selfdemo; return; end
if nargin<5, text4record=''; end
if nargin<6, forcedRecording=0; end

if forcedRecording
	message=sprintf('\tPress "Enter" to record "%s" for %g-seconds: ', text4record, duration);
	if isempty(text4record), message=sprintf('\tPress "Enter" to record for %g-sec: ', duration); end
	fprintf('%s', message);
	pause
	userInput = [];
	fprintf('\n');
else
	message=sprintf('\tPress "Enter" and record "%s" for %g-seconds (Or press "n" before "Enter" to skip): ', text4record, duration);
	if isempty(text4record), message=sprintf('\tPress "Enter" and record for %g-sec: ', duration); end
%	fprintf('%s', message);
	userInput = input(message, 's');		% 分開成兩部分，才能完整顯示特殊中文，如「淚」、「許」等。
end

if isempty(userInput),
	fprintf('\tRecording... ');
	if verLessThan('matlab', '8.1')
		if nbits==8
			y=wavrecord(duration*fs, fs, 'uint8');
		elseif nbits==16
			y=wavrecord(duration*fs, fs, 'int16');
		else
			error('Wrong nbits');
		end
		y = double(y);			% Convert from a uint8 to double array
		y = y/(2^nbits/2);		% Make y zero mean and unity maximum
	else
		recorder = audiorecorder(fs, nbits, 1);
		recordblocking(recorder, duration);
		y = getaudiodata(recorder);	% "double" by default
	end
	y = y-mean(y);
	plot((1:length(y))/fs, y);
	waveFile=strrep(waveFile, '\', '/');
	title(strPurify(['Waveform of "', waveFile, '"']));
	axis([-inf inf -1 1]);
	figure(gcf);
	if verLessThan('matlab', '8.1')
		wavwrite(y, fs, nbits, waveFile);
	else
		audiowrite(waveFile, y, fs, 'BitsPerSample', nbits);
	end
	fprintf('Finish recording and save to "%s".\n', waveFile);
	fprintf('\tPlayback...\n'); sound(y, fs);
	success=1;
else
	fprintf('\tNothing is saved.\n');
	success=0;
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
