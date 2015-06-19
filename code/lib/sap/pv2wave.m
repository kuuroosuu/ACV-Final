function y=pv2wave(pv, pitchRate, fs, method, plotOpt)
% pv2wave: Pitch vector (pv) to wave signals conversion
%	
%	Usage: pv2wave(pv, pitchRate, fs, method, plotOpt)
%
%	Example:
%		pv=[60.3304 59.2131 58.6804 59.2131 59.7627 60.3304 60.3304 60.3304 59.7627 59.7627 59.7627 59.7627 0 0 0 60.3304 59.2131 58.6804 58.6804 59.2131 59.2131 59.7627 59.7627 59.7627 59.2131 59.7627 59.7627 61.5248 63.4868 65.6999 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 67.35 66.5053 0 0 0 0 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 67.35 67.35 67.35 67.35 68.238 69.174 69.174 68.238 68.238 68.238 68.238 68.238 68.238 68.238 68.238 0 0 0 0 69.174 68.238 68.238 68.238 68.238 68.238 68.238 69.174 68.238 68.238 68.238 0 0 0 0 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 66.5053 67.35 67.35 65.699915 63.4868 62.8078 62.8078 63.4868 64.9304 64.1935 64.1935 64.1935 64.9304 64.9304 64.9304 64.1935 64.9304 64.9304 64.9304 0 0 0 0 64.9304 64.1935 64.1935 64.1935 64.1935 64.9304 64.9304 64.9304 64.9304 64.9304 64.1935 64.1935 63.486821 64.1935 63.4868 63.4868 63.4868 63.4868 63.4868 63.4868 63.4868 63.4868 64.1935 63.4868 0 0 63.4868 64.1935 63.4868 63.4868 63.4868 63.4868 63.4868 63.4868 63.4868 64.1935 63.4868 63.4868 0 0 0 0 60.9173 60.9173 60.9173 61.5248 61.5248 61.5248 61.5248 61.5248 61.5248 61.5248 62.8078 0 0 61.5248 60.9173 60.9173 60.9173 61.5248 61.5248 61.5248 61.5248 61.5248 61.5248 61.5248 62.1544 0 0 0 0 59.7627 59.2131 59.2131 59.7627 59.7627 59.7627 59.7627 59.2131 59.2131 59.2131 59.2131 59.2131 59.7627 60.3304 60.3304];
%		fs=8000;
%		pitchRate=8000/256;
%		method=1;
%		plotOpt=1;
%		y=pv2wave(pv, pitchRate, fs, method, plotOpt);
%		sound(y, fs);
%
%	See also pvPlay, pv2waveMex.

%	Roger Jang, 20041019, 20070531

if nargin<1, selfdemo; return; end
if nargin<2, pitchRate=8000/256; end
if nargin<3, fs=8000; end
if nargin<4, method=1; end
if nargin<5, plotOpt=0; end

switch (method)
	case 1
		pv=pv(:)';
		pv = double(pv);
		pv(isnan(pv)) = 0;    % Convert NaN to 0
		pvLen=length(pv);
		duration=ones(1, pvLen)/pitchRate;
	%	[note(1:length(pvLen)).pitch] = deal(pv);
	%	[note(1:length(pvLen)).duration] = deal(duration);
		note=[pv; duration]; note=note(:);
		y=note2wave(note, 1, fs, 2, plotOpt);
	case 2
		index=find(pv==0);
		freq=440*2.^((pv-69)/12);	% pitch=69+12*log2(freq/440)
		freq(index)=0;
		frameSize=round(fs/pitchRate);
		y=zeros(1, frameSize*length(pv));
		for i=1:length(freq)
			t=(1:frameSize)/fs;
			if freq(i)==0
				frame=0*t;
			else
				frame=sin(2*pi*freq(i)*t).*((frameSize-1:-1:0)/(frameSize-1));
			end
			y((i-1)*frameSize+1:i*frameSize)=frame;
		end
	otherwise
		fprintf('Unknown method in pv2wave()!\n');
end

if plotOpt
	temp=pv;
	temp(temp==0)=nan;
	subplot(2,1,1);
	plot((1:length(temp))/pitchRate, temp, (1:length(temp))/pitchRate, temp, '.');
	ylabel('Pitch (semitone)');
	subplot(2,1,2);
	plot((1:length(y))/fs, y);
	ylabel('Waveform');
	xlabel('Time (sec)');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
