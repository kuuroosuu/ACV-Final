function outWave=beepCreate(beatTime, fs, origWave)
%beepCreate: Create a sequence of beeps indicating the location of beats, etc.
%	Usage:
%		outWave=beepCreate(beatTime, fs)
%		outWave=beepCreate(beatTime, fs, origWave)
%
%	Description:
%		outWave=beepCreate(beatTime, fs) returns a sequence of beeps indicating the location of beats.
%		outWave=beepCreate(beatTime, fs, origWave) add the beeps to the origWave.

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, origWave=zeros(max(beatTime)*fs+2000, 1); end

% Create a short beep with duration 100ms and frequency 2000Hz
beepDuration=0.2;
beepFreq=2000;
time=(0:round(beepDuration*fs))';
beep=time.*exp(-time/((beepDuration*fs)/10)).*cos(2*pi*time/fs*beepFreq);
beep=beep/max(abs(beep));	% Normalization

outWave=origWave;
beatSample=round(beatTime*fs);
for i=1:length(beatTime)
	index1=round(beatSample(i));
	index2=min(index1+length(beep)-1, length(outWave));
	outWave(index1:index2)=outWave(index1:index2)+beep(1:index2-index1+1);
end

% ====== Self demo
function selfdemo
beatTime=[0.1, 0.2, 0.4, 0.5, 0.6];
fs=16000;
beepTrain=beepCreate(beatTime, fs);
sound(beepTrain, fs);
plot((1:length(beepTrain))/fs, beepTrain);
