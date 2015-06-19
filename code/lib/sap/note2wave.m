function wave=note2wave(note, timeUnit, fs, method, plotOpt)
%note2wave: Music note to wave conversion
%	Usage: wave=note2wave(note, timeUnit, fs, method, plotOpt)
%		note: music note vector, with two possible formats:
%			Format 1:
%				note.pitch: pitch in semitones, where zero/nan is considered as silence
%				note.duration: duration in timeUnit
%			Format 2:
%				note = [pitch1, duration1, pitch2, duration2, ...]
%		fs: sample rate of the output wave signals
%		method: method for generating the wave signals
%			method=1 ===> Piano-like energy profile
%			method=2 ===> Flat energy profile
%		plotOpt: 1 for plotting
%		wave: the output wave signals
%
%	Example:
%		pitch = [-2 -5 -5 -4 -7 -7 -9 -7 -5 -4 -2 -2 -2]+69;
%		duration = [1   1  2  1  1  2  1  1  1  1  1  1  2]/2;
%		note=[pitch; duration]; note=note(:)';
%		timeUnit=1;
%		fs=16000;
%		fprintf('Using method 1 in note2wave.m. Hit return to start...'); pause; fprintf('\n');
%		method=1; wave=note2wave(note, timeUnit, fs, method); sound(wave, fs');
%		fprintf('Using method 2 in note2wave.m. Hit return to start...'); pause; fprintf('\n');
%		method=2; wave=note2wave(note, timeUnit, fs, method); sound(wave, fs);

%	Roger Jang, 19990516, 20041020, 20060303, 20070531

if nargin<1, selfdemo; return; end
if nargin<2, timeUnit=1; end
if nargin<3, fs=16000; end
if nargin<4, method=1; end
if nargin<5, plotOpt=0; end

if isstruct(note)
	pitch=[note.pitch];
	duration=[note.duration]*timeUnit;
else
	pitch=note(1:2:end);
	duration=note(2:2:end)*timeUnit;
end
pitch(isnan(pitch))=0;

switch (method)
	case 1
		freq = 440*2.^((pitch-69)/12);
		wave = [];
		opt=noteEnvelope('defaultOpt');
		for i=1:length(pitch)
			t = 0:1/fs:duration(i);		% In seconds
			if pitch(i)~=0
				y = sin(2*pi*freq(i)*t).*noteEnvelope(t, opt);
			%	y = sin(2*pi*freq(i)*t);
			else
				y = 0*t;
			end
			wave = [wave(1:end-1) y];          
		end
	case 2
		wave=note2waveMex(pitch, duration, fs);
	otherwise
		fprintf('Unknown method in note2wave.m!\n');
end

if plotOpt
	subplot(2,1,1);
	notePlot(note);
	subplot(2,1,2);
	plot((1:length(wave))/fs, wave);
	xlabel('Time (Seconds)'); ylabel('Waveform');
	set(gca, 'xlim', [-inf, inf]);
	drawnow;
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
